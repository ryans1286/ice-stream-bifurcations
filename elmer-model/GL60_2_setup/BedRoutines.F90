MODULE BedRoutines

  USE Types
  USE SParIterComm
  USE MainUtils

  IMPLICIT NONE

  CONTAINS

 FUNCTION PointInPolygon(Polygon, Point, buffer, ZeroFlag) RESULT(inside)
    REAL(kind=dp) :: polygon(:,:)
    REAL(kind=dp) :: point(2)
    REAL(kind=dp), OPTIONAL :: buffer
    LOGICAL, OPTIONAL :: ZeroFlag
    !--------------------------------
    REAL(kind=dp), ALLOCATABLE :: ZPolygon(:,:)
    INTEGER :: n,i,windingnumber
    REAL(kind=dp) :: left,bounding(4),dist,ZPoint(2)
    LOGICAL :: inside,Zero,OutOfBounds

    IF(SIZE(polygon(:,1)) /= 2) CALL FATAL('PointInPolygon2D', 'Please provide a 2D array with x and y coords')

    IF(PRESENT(ZeroFlag)) THEN
      Zero = ZeroFlag
    ELSE
      Zero = .TRUE.
    END IF
    n=SIZE(polygon(1,:))

    IF(Zero) THEN
      ZPoint = Point
      ALLOCATE(ZPolygon(2,n))
      ZPolygon = Polygon
      CALL ZeroPolygon(Polygon, Point)
    END IF

    windingnumber=100
    DO i=1, n-1
      ! polygon y i <= point y
      IF(Polygon(2,i) <= Point(2)) THEN !start with y<=P.y
        IF(Polygon(2, i+1) > Point(2)) THEN !upward crossing
          left=IsLeft(Polygon(:, i), Polygon(:, i+1), Point(:))
          IF(left > 0.0_dp) THEN !p is to left of intersect
            windingnumber=windingnumber+1 !valid up intersect
          END IF
        END IF
      ELSE    !start at y> point y
        IF(Polygon(2, i+1) <= Point(2)) THEN ! downward crossing
          Left = IsLeft(Polygon(:, i), Polygon(:, i+1), Point(:))
          IF(left < 0.0_dp) THEN ! p right of edge
            windingnumber=windingnumber-1
          END IF
        END IF
      END IF
    END DO

    IF(windingnumber /= 100) THEN
      inside = .TRUE.
    ELSE
      inside = .FALSE.
    END IF

    IF(.NOT. PRESENT(buffer) .OR. inside) THEN
      IF(Zero) THEN
        Polygon = ZPolygon
        Point = ZPoint
      END IF
      RETURN
    END IF


    bounding(1) = MINVAL(Polygon(1,:)) - buffer
    bounding(2) = MAXVAL(Polygon(1,:)) + buffer
    bounding(3) = MINVAL(Polygon(2,:)) - buffer
    bounding(4) = MAXVAL(Polygon(2,:)) + buffer

    OutOfBounds = .FALSE.
    IF(Point(1) <  bounding(1)) THEN
      OutOfBounds = .TRUE.
    ELSE IF(Point(1) > bounding(2)) THEN
      OutOfBounds = .TRUE.
    ELSE IF(Point(2) < bounding(3)) THEN
      OutOfBounds = .TRUE.
    ELSE IF(Point(2) > bounding(4)) THEN
      OutOfBounds = .TRUE.
    END IF

    IF(OutOfBounds) THEN
      IF(Zero) THEN
        Polygon = ZPolygon
        Point = ZPoint
      END IF
      RETURN
    END IF

    DO i=1, n-1
      dist = PointLineSegmDist2D(Polygon(:,i), Polygon(:,i+1),Point)
      IF(dist < buffer) THEN
        windingnumber = 99
        EXIT
      END IF
    END DO

    IF(windingnumber /= 100) THEN
      inside = .TRUE.
    ELSE
      inside = .FALSE.
    END IF

    IF(Zero) THEN
      Polygon = ZPolygon
      Point = ZPoint
    END IF

  END FUNCTION PointInPolygon

    SUBROUTINE ZeroPolygon(Polygon, Point)
      REAL(kind=dp) :: Polygon(:,:), Point(2)
      REAL(kind=dp) :: minx, miny

      minx = MINVAL(Polygon(1,:))
      miny = MINVAL(Polygon(2,:))

      Polygon(1,:) = Polygon(1,:) - minx
      Polygon(2,:) = Polygon(2,:) - miny

      Point(1) = Point(1) - minx
      Point(2) = Point(2) - miny

    END SUBROUTINE ZeroPolygon

    FUNCTION IsLeft(a, b, c) RESULT(d)
      REAL(kind=dp) :: a(2), b(2), c(2), d

      d = (b(1)-a(1)) * (c(2)-a(2)) - &
      (c(1)-a(1)) * (b(2)-a(2))

    END FUNCTION Isleft

    FUNCTION  PointLineSegmDist2D(a, b, c)  RESULT (pdis)
      REAL(KIND=dp) :: a(2), b(2), c(2), n(2), v(2), dd, t, pdis
      n=b-a                      ! Vector ab
      dd = (n(1)**2.+n(2)**2.)   ! Length of ab squared
      dd = DOT_PRODUCT(n,n) ! alternative writing
      t = DOT_PRODUCT(c-a,b-a)/dd
      dd = MAXVAL( (/0.0_dp, MINVAL( (/1.0_dp,t/) ) /) )
      v = c - a - dd * n
      pdis=sqrt(v(1)**2.+v(2)**2.)
    END FUNCTION PointLineSegmDist2D


   SUBROUTINE GetFrontCornersCopy(Model, Solver, FrontMaskName, LeftMaskName, RightMaskName, FrontLeft, FrontRight)

    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    CHARACTER(LEN=*):: FrontMaskName, LeftMaskName, RightMaskName
    !--------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: NullSolver => NULL(), AdvSolver
    TYPE(Valuelist_t), POINTER :: SolverParams, AdvParams
    INTEGER :: i,j,k, dummyint, LeftRoot, RightRoot, ierr, NNodes,RCounter, LCounter,&
        Nl,Nr, Naux, ok, RightTotal, LeftTotal, Nrail, CornersTotal, Counter, side
    REAL(KIND=dp) :: FrontLeft(2), FrontRight(2), buffer, xx, yy, mindist, tempdist
    INTEGER, POINTER :: FrontPerm(:)=>NULL(), TopPerm(:)=>NULL(), &
        LeftPerm(:)=>NULL(), RightPerm(:)=>NULL(), SidePerm(:)
    LOGICAL :: FoundRight, FoundLeft, reducecorners(2), Found
    LOGICAL, ALLOCATABLE :: PFoundRight(:), PFoundLeft(:), InFront(:), Duplicate(:)
    INTEGER, ALLOCATABLE ::  PRightCount(:), PLeftCount(:), disps(:),&
        PCount(:), jmin(:), Corner(:)
    REAL(kind=dp), ALLOCATABLE :: xL(:),yL(:),xR(:),yR(:), xRail(:), yRail(:),&
        AllCorners(:), PAllCorners(:), MinDists(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: TopMaskName,SolverName = "GetFrontCorners",&
         RightRailFName, LeftRailFName, Adv_EqName
    INTEGER, PARAMETER :: io=20

    NNodes = Model % Mesh % NumberOfNodes
    Mesh => Model % Mesh
    SolverParams => Solver % Values

    ALLOCATE(FrontPerm(NNodes), TopPerm(NNodes), LeftPerm(NNodes),&
        RightPerm(NNodes))
    !FrontMaskName = "Random Out Mask"
    TopMaskName = "Top Surface Mask"
    CALL MakePermUsingMask( Model, Solver, Mesh, TopMaskName, &
      .FALSE., TopPerm, dummyint)
    CALL MakePermUsingMask( Model, Solver, Mesh, FrontMaskName, &
      .FALSE., FrontPerm, dummyint)
    !LeftMaskName = "Random Left Sidewall Mask"
    !RightMaskName = "Random Right Sidewall Mask"
    !Generate perms to quickly get nodes on each boundary
    CALL MakePermUsingMask( Model, Solver, Mesh, LeftMaskName, &
      .FALSE., LeftPerm, dummyint)
    CALL MakePermUsingMask( Model, Solver, Mesh, RightMaskName, &
         .FALSE., RightPerm, dummyint)

    PRINT*, COUNT(LeftPerm>0), COUNT(RightPerm > 0), COUNT(FrontPErm>0), COUNT(TopPerm > 0)

    FoundLeft=.FALSE.
    FoundRight=.FALSE.
    RCounter= 0; LCounter=0
    DO i=1,NNodes
       IF( (TopPerm(i) >0 ) .AND. (FrontPerm(i) >0 )) THEN
         IF( LeftPerm(i) >0  ) THEN
            FrontLeft(1) = Mesh % Nodes % x(i)
            FrontLeft(2) = Mesh % Nodes % y(i)
            LCounter = LCounter + 1
            FoundLeft = .TRUE.
         ELSE IF ( RightPerm(i) >0  ) THEN
            FrontRight(1) = Mesh % Nodes % x(i)
            FrontRight(2) = Mesh % Nodes % y(i)
            RCounter = RCounter + 1
            FoundRight = .TRUE.
         END IF
       END IF
    END DO

    PRINT*, ParENv % MYPE, FoundRight, FoundLeft

    ALLOCATE(PFoundRight(ParEnv % PEs), PFoundLeft(ParEnv % PEs))
    CALL MPI_ALLGATHER(FoundRight, 1, MPI_LOGICAL, PFoundRight, 1, &
            MPI_LOGICAL, ELMER_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER(FoundLeft, 1, MPI_LOGICAL, PFoundLeft, 1, &
            MPI_LOGICAL, ELMER_COMM_WORLD, ierr)

    DO i=1, ParEnv % PEs
      IF(.NOT. PFoundLeft(i) .AND. .NOT. PFoundRight(i)) CYCLE
      IF(PFoundLeft(i)) LeftRoot = i-1
      IF(PFoundRight(i)) RightRoot = i-1
    END DO

    IF(ALL(.NOT. PFoundLeft)) CALL FATAL(SolverName, 'Unable to find left corner')
    IF(ALL(.NOT. PFoundRight)) CALL FATAL(SolverName, 'Unable to find right corner')

    CALL MPI_BCAST(FrontLeft,2,MPI_DOUBLE_PRECISION, LeftRoot, ELMER_COMM_WORLD, ierr)

    CALL MPI_BCAST(FrontRight,2,MPI_DOUBLE_PRECISION, RightRoot, ELMER_COMM_WORLD, ierr)
    DEALLOCATE(FrontPerm, TopPerm, LeftPerm, RightPerm)

  END SUBROUTINE GetFrontCornersCopy

  FUNCTION GetFrontOrientationCopy(Model, FrontMaskName, LeftMaskName, RightMaskName) RESULT (Orientation)
    TYPE(Model_t) :: Model
    TYPE(Mesh_t),POINTER :: Mesh
    CHARACTER(LEN=*):: FrontMaskName, LeftMaskName, RightMaskName
    !--------------------------
    TYPE(Solver_t), POINTER :: NullSolver => NULL()
    TYPE(Variable_t), POINTER :: TimeVar
    INTEGER :: i,dummyint,FaceNodeCount, ierr, proc
    REAL(KIND=dp) :: Orientation(3),OrientSaved(3),xLeft,yLeft,xRight,yRight
    REAL(KIND=dp) :: RecvXL,RecvYL,RecvXR,RecvYR,Temp,PrevTime
    REAL(KIND=dp), POINTER :: PArray(:,:) => NULL()
    INTEGER, POINTER :: Perm(:), FrontPerm(:)=>NULL(), TopPerm(:)=>NULL(), &
        FrontNodeNums(:)=>NULL(),LeftPerm(:)=>NULL(), RightPerm(:)=>NULL()
    LOGICAL :: FirstTime=.TRUE.,Constant,Debug=.TRUE.,Parallel,&
         HaveRight,HaveLeft, Boss, FirstThisTime
    CHARACTER(LEN=MAX_NAME_LEN) :: TopMaskName
    INTEGER :: status(MPI_STATUS_SIZE), iLeft, iRight
    SAVE :: FirstTime,Constant,PArray,Parallel, Boss, FirstThisTime
    SAVE :: PrevTime

    Parallel = (ParEnv % PEs > 1)
    Boss = (ParEnv % MyPE == 0) .OR. (.NOT. Parallel)
    Constant = .FALSE.
    FirstThisTime=.TRUE.

    ! check whether already did a front orientation computation this timestep
    ! Changed Model % Mesh % Variables to avoid segfault as when calling vtusolver after mmg step as
    ! Model % Variables lost after vtuoutput 
    IF(Constant .OR. (.NOT. FirstThisTime) ) THEN
      Orientation = OrientSaved
      RETURN
    ELSE
      PRINT *, 'computing orientation'
      Orientation(3) = 0.0_dp ! always set z-component to 0
      Mesh => Model % Mesh
      !Get the front line
      !FrontMaskName = "Calving Front Mask"
      TopMaskName = "Top Surface Mask"
      CALL MakePermUsingMask( Model, NullSolver, Mesh, TopMaskName, &
        .FALSE., TopPerm, dummyint)
      CALL MakePermUsingMask( Model, NullSolver, Mesh, FrontMaskName, &
        .FALSE., FrontPerm, FaceNodeCount)
      !LeftMaskName = "Left Sidewall Mask"
      !RightMaskName = "Right Sidewall Mask"
      !Generate perms to quickly get nodes on each boundary
      CALL MakePermUsingMask( Model, NullSolver, Mesh, LeftMaskName, &
        .FALSE., LeftPerm, dummyint)
      CALL MakePermUsingMask( Model, NullSolver, Mesh, RightMaskName, &
           .FALSE., RightPerm, dummyint)
      iLeft=0
      iRight=0
      HaveLeft=.FALSE.
      HaveRight=.FALSE.
      DO i=1,Mesh % NumberOfNodes
         IF( (TopPerm(i) >0 ) .AND. (FrontPerm(i) >0 )) THEN
           IF( LeftPerm(i) >0  ) THEN
              xLeft = Mesh % Nodes % x(i)
              yLeft = Mesh % Nodes % y(i)
              HaveLeft =.TRUE.
           ELSE IF ( RightPerm(i) >0  ) THEN
              xRight = Mesh % Nodes % x(i)
              yRight = Mesh % Nodes % y(i)
              HaveRight =.TRUE.
           END IF
         END IF
      END DO
      IF (Debug)  PRINT *, 'GetFrontOrientation: HaveLeft, HaveRight', HaveLeft, HaveRight
      IF (Parallel) THEN
         IF (HaveLeft) PRINT *, 'GetFrontOrientation: xL, yL', xLeft, yLeft
         IF (HaveRight)  PRINT *, 'GetFrontOrientation: xR, yR', xRight, yRight
         IF (Debug) PRINT *, 'communicate the corners'
         IF (HaveLeft  .AND. (ParEnv % MyPE>0)) THEN ! left not in root
            iLeft=ParEnv % MyPE
            CALL MPI_BSEND(xLeft, 1, MPI_DOUBLE_PRECISION, &
                 0 ,7001, ELMER_COMM_WORLD, ierr )
            CALL MPI_BSEND(yLeft, 1, MPI_DOUBLE_PRECISION, &
                 0 ,7002, ELMER_COMM_WORLD, ierr )
            IF (Debug) PRINT *, 'sending left'
         END IF
         IF (HaveRight .AND. (ParEnv % MyPE>0) ) THEN ! right not in root
            iRight=ParEnv % MyPE
            CALL MPI_BSEND(xRight, 1, MPI_DOUBLE_PRECISION, &
                 0 , 7003, ELMER_COMM_WORLD, ierr )
            CALL MPI_BSEND(yRight, 1, MPI_DOUBLE_PRECISION, &
                 0 , 7004, ELMER_COMM_WORLD, ierr )
            IF (Debug) PRINT *, 'sending right'
         END IF
         IF (Debug) PRINT *, 'sent the corners'
         IF (Boss) THEN
            IF (Debug) PRINT *, ParEnv % PEs
            IF (.NOT.HaveLeft) THEN
                  CALL MPI_RECV(RecvXL,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
                       7001,ELMER_COMM_WORLD, status, ierr )
                  CALL MPI_RECV(RecvYL,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
                       7002,ELMER_COMM_WORLD, status, ierr )
                  xLeft=RecvXL
                  yLeft=RecvYL
            END IF
            IF (.NOT. HaveRight) THEN
                  CALL MPI_RECV(RecvXR,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
                       7003,ELMER_COMM_WORLD, status, ierr )
                  CALL MPI_RECV(RecvYR,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
                       7004,ELMER_COMM_WORLD, status, ierr )
                  xRight=RecvXR
                  yRight=RecvYR
            END IF
            IF (Debug) PRINT *, 'received corners'
            IF (Debug) PRINT *, 'GetFrontOrientation: Boss xL, yL, xR, yR', xLeft, yLeft, xRight, yRight
         END IF
      END IF ! end if parallel
      IF (Boss) THEN ! root or not parallel
      IF( ABS(xLeft-xRight) < AEPS) THEN
         ! front orientation is aligned with y-axis
         Orientation(2) =  0.0_dp
         IF(yRight > yLeft) THEN
            Orientation(1)=1.0_dp
         ELSE
            Orientation(1)=-1.0_dp
         END IF
      ELSE IF (ABS(yLeft-yRight)<AEPS) THEN
         ! front orientation is aligned with x-axis
         Orientation(1) = 0.0_dp
         IF(xRight > xLeft) THEN
            Orientation(2)=1.0_dp
         ELSE
            Orientation(2)=-1.0_dp
         END IF
      ELSE
         ! set dot product equal to 0
         ! no need to ensure it is unit normal, done in ComputeRotation
         IF(xRight > xLeft) THEN
            Orientation(2)=1.0_dp
         ELSE
            Orientation(2)=-1.0_dp
         END IF
         Orientation(1)=Orientation(2)*(yRight-yLeft)/(xLeft-xRight)
      END IF
      END IF !boss
      IF (Parallel) CALL MPI_BCAST(Orientation,3,MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
      ! deallocations
       DEALLOCATE(FrontPerm, TopPerm, LeftPerm, RightPerm)
    END IF
    Temp=(Orientation(1)**2+Orientation(2)**2+Orientation(3)**2)**0.5
    Orientation=Orientation/Temp ! normalized the orientation
    IF((.NOT. Constant).AND.Debug)  PRINT *, "GetFrontOrientation: ", Orientation,'part',ParEnv % MyPE
    FirstThisTime=.FALSE.
  END FUNCTION GetFrontOrientationCopy

END MODULE BedRoutines

SUBROUTINE CreateBed( Model, Solver, dt, TransientSimulation )
  USE CalvingGeometry
  USE MainUtils
  USE BedRoutines
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!-------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: NullSolver => NULL()
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  INTEGER :: i,NNodes,dummyint,ierr
  INTEGER, POINTER :: Perm(:), RightPerm(:)=>NULL(),LeftPerm(:)=>NULL()
  REAL(KIND=dp), POINTER :: Values(:)
  REAL(KIND=dp) :: LeftOrientation(3),LRotationMatrix(3,3),LUnRotationMatrix(3,3),&
    RightOrientation(3),RRotationMatrix(3,3),RUnRotationMatrix(3,3),Slope1,Slope2,&
    elevation,x,y,z0,width,rx,ry,NodeHolder(3),rz0,FrontLL(2),FrontLR(2),&
    FrontRL(2),FrontRR(2),RotFrontLL(2),RotFrontLR(2),RotFrontRL(2),RotFrontRR(2),rightnodewidth,&
    addition,leftnodewidth,ydif,extra_slope,extra_z,slopedist,zdif,flatdistance,dist_rx,&
    left_lengthtobif,right_lengthtobif,left_yorigin,right_yorigin,cross(2),LeftPoly(2,5),RightPoly(2,5),&
    m1,m2,left_right_min,left_left_min,right_left_min,right_right_min
  LOGICAL :: Found,inside
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,LeftLeftName,LeftRightName,LeftFrontName,RightLeftName,&
    RightRightName,RightFrontName

  SolverName = "CreateBed"

  Mesh => Model % Mesh
  NNodes = Mesh % NumberOfNodes

  Values => Solver % Variable % Values
  Perm => Solver % Variable % Perm

  Params => Solver % Values

  ! read in variables
  Slope1 = ListGetConstReal(Params, "Main trunk slope",Found, UnfoundFatal=.TRUE.)
  Slope2 = ListGetConstReal(Params, "Second trunk slope",Found, UnfoundFatal=.TRUE.)
  Width = ListGetConstReal(Params, "Main trunk width",Found, UnfoundFatal=.TRUE.)
  z0 = ListGetConstReal(Params, "z0",Found, UnfoundFatal=.TRUE.)
  flatdistance = ListGetConstReal(Params, "Second trunk flat distance",Found, UnfoundFatal=.TRUE.)

  LeftLeftName = ListGetString(Params,"Left Channel Left Name", UnfoundFatal=.TRUE.)
  LeftRightName = ListGetString(Params,"Left Channel Right Name", UnfoundFatal=.TRUE.)
  LeftFrontName = ListGetString(Params,"Left Channel Front Name", UnfoundFatal=.TRUE.)

  RightLeftName = ListGetString(Params,"Right Channel Left Name", UnfoundFatal=.TRUE.)
  RightRightName = ListGetString(Params,"Right Channel Right Name", UnfoundFatal=.TRUE.)
  RightFrontName = ListGetString(Params,"Right Channel Front Name", UnfoundFatal=.TRUE.)

  ! get left orientation
  LeftOrientation = GetFrontOrientationCopy(Model, LeftFrontName, LeftLeftName, LeftRightName)
  LRotationMatrix = ComputeRotationMatrix(LeftOrientation)
  LUnRotationMatrix = TRANSPOSE(LRotationMatrix)

  ! get right orientation
  RightOrientation = GetFrontOrientationCopy(Model, RightFrontName, RightLeftName, RightRightName)
  RRotationMatrix = ComputeRotationMatrix(RightOrientation)
  RUnRotationMatrix = TRANSPOSE(RRotationMatrix)

  ! get front corners
  CALL GetFrontCornersCopy(Model, Solver, LeftFrontName, LeftLeftName, LeftRightName, FrontLL, FrontLR)
  CALL GetFrontCornersCopy(Model, Solver, RightFrontName, RightLeftName, RightRightName, FrontRL, FrontRR)

  ! orientate front corners
  NodeHolder(1:2) = FrontLR; NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(LRotationMatrix,NodeHolder)
  RotFrontLR = NodeHolder(2:3)

  NodeHolder(1:2) = FrontLL; NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(LRotationMatrix,NodeHolder)
  RotFrontLL = NodeHolder(2:3)
  leftnodewidth = ABS(RotFrontLL(1)-RotFrontLR(1))

  NodeHolder(1:2) = FrontRR; NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(RRotationMatrix,NodeHolder)
  RotFrontRR = NodeHolder(2:3)

  NodeHolder(1:2) = FrontRL; NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(RRotationMatrix,NodeHolder)
  RotFrontRL = NodeHolder(2:3)
  rightnodewidth = ABS(RotFrontRL(1)-RotFrontRR(1))

  !get min point of z origin in rotated coords

  !left first
  CALL MakePermUsingMask( Model, Solver, Mesh, LeftRightName, &
    .FALSE., RightPerm, dummyint)
  CALL MakePermUsingMask( Model, Solver, Mesh, LeftLeftName, &
    .FALSE., LeftPerm, dummyint)

  left_lengthtobif = HUGE(1.0_dp);
  left_left_min = HUGE(1.0_dp); left_right_min = HUGE(1.0_dp)
  DO i=1,NNodes
    IF(RightPerm(i) == 0 .AND. LeftPerm(i) == 0) CYCLE
    x = Mesh % Nodes % x(i)
    y = Mesh % Nodes % y(i)

    NodeHolder(1) = Mesh % Nodes % x(i)
    NodeHolder(2) = Mesh % Nodes % y(i)
    NodeHolder(3) = 0.0_dp

    NodeHolder = MATMUL(LRotationMatrix,NodeHolder)
    ry = NodeHolder(3)

    IF(LeftPerm(i) > 0) THEN
      IF(ry < left_left_min) left_left_min = ry
    END IF

    IF(RightPerm(i) > 0) THEN
      IF(ry < left_right_min) left_right_min = ry
    END IF

    IF(LeftPerm(i) == 0) CYCLE
    IF(y < left_lengthtobif) left_lengthtobif = y
  END DO

  ! parallel comm
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,left_lengthtobif, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,left_right_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,left_left_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)

  left_yorigin = left_right_min!MIN(left_left_min, left_right_min)

  PRINT*, 'mins', left_left_min, left_right_min, 'origin', left_yorigin, 'dif', left_left_min-left_right_min

  DEALLOCATE(LeftPerm,RightPerm)

  !now right
  CALL MakePermUsingMask( Model, Solver, Mesh, RightRightName, &
    .FALSE., RightPerm, dummyint)
  CALL MakePermUsingMask( Model, Solver, Mesh, RightLeftName, &
    .FALSE., LeftPerm, dummyint)

  right_lengthtobif = HUGE(1.0_dp); cross = HUGE(1.0_dp)
  right_left_min = HUGE(1.0_dp); right_right_min = HUGE(1.0_dp)
  DO i=1,NNodes
    IF(RightPerm(i) == 0 .AND. LeftPerm(i) == 0) CYCLE
    x = Mesh % Nodes % x(i)
    y = Mesh % Nodes % y(i)

    NodeHolder(1) = Mesh % Nodes % x(i)
    NodeHolder(2) = Mesh % Nodes % y(i)
    NodeHolder(3) = 0.0_dp

    NodeHolder = MATMUL(RRotationMatrix,NodeHolder)
    ry = NodeHolder(3)

    IF(LeftPerm(i) > 0) THEN
      IF(ry < right_left_min) right_left_min = ry
    END IF

    IF(RightPerm(i) > 0) THEN
      IF(ry < right_right_min) right_right_min = ry
    END IF

    IF(x<cross(1) .AND. LeftPerm(i) > 0) THEN
      cross(1) = x
      cross(2) = y
    END IF

    IF(RightPerm(i) == 0) CYCLE
    IF(y < right_lengthtobif) right_lengthtobif = y
  END DO

  ! parallel comm
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,right_lengthtobif, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,cross, 2, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,right_right_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,right_left_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)

  right_yorigin = right_left_min!MIN(right_left_min, right_right_min)

  DEALLOCATE(LeftPerm,RightPerm)

  IF(ABS(left_lengthtobif-right_lengthtobif) > AEPS) CALL FATAL(SolverName, 'Bifurcation start at different trunk lengths')

  ! create polygon of left birfucation
  LeftPoly(1,1) = 0.0_dp; LeftPoly(2,1) = left_lengthtobif
  LeftPoly(:,2) = FrontLL
  LeftPoly(:,3) = FrontLR
  LeftPoly(:,4) = cross
  LeftPoly(:,5) = LeftPoly(:,1)

  ! create polygon of left birfucation
  RightPoly(1,1) = width; RightPoly(2,1) = right_lengthtobif
  RightPoly(:,2) = FrontRR
  RightPoly(:,3) = FrontRL
  RightPoly(:,4) = cross
  RightPoly(:,5) = RightPoly(:,1)

  IF(ParEnv % MYPE == 0 .OR. ParEnv % MYPE == 2) THEN
    PRINT*, ParEnv % MYPE, 'leftpoly', LeftPoly
    PRINT*, ParEnv % MYPE, 'rightpoly', RightPoly
    PRINT*, ParEnv % MyPE, 'cross', cross
  END IF

  ! loop through nodes assigning values
  DO i=1, NNodes
    x = Mesh % Nodes % x(i)
    y = Mesh % Nodes % y(i)
    
    IF(y <= left_lengthtobif) THEN ! on main trunk
      elevation = (z0 - slope1*y)
      extra_z=0.0_dp
      zdif=0.0_dp
    ELSE ! on bifurications or join

      IF(x < cross(1)) THEN ! left

        inside = PointInPolygon(LeftPoly, (/x,y/), 0.1_dp, .TRUE.)

        IF(inside) THEN

!#if 0
          !rotatated coord system with origin as start of second trunk
          NodeHolder(1) = Mesh % Nodes % x(i)
          NodeHolder(2) = Mesh % Nodes % y(i)
          NodeHolder(3) = 0.0_dp

          NodeHolder = MATMUL(LRotationMatrix,NodeHolder)

          rx = leftnodewidth - (NodeHolder(2) - RotFrontLL(1))
          ry = NodeHolder(3) - left_yorigin

          ! distance along divide between trunks
          !addition = (rx**2+(rx*LeftOrientation(2))**2)**0.5
          !m = (cross(2)-left_lengthtobif)/(cross(1)-width)
          !addition = -1*m1*(leftnodewidth-rx)

          ! distance in flow direction since divide
          !ydif = ry - rx*LeftOrientation(2)
          m1 = (cross(2)-left_lengthtobif) / leftnodewidth
          m2 = (left_left_min-left_right_min) / leftnodewidth
          ydif =  ry - rx * m2 !- left_y
          addition = m1*(leftnodewidth-rx)


          ! height at divide
          rz0 = (z0 - slope1 * (left_lengthtobif+addition))

          ! additional slope twisting. Zdiff is the additional heightloss required to flatten out
          ! second trunk. It is distance above lowest point of ice divide minus height loss from extra length
          dist_rx = rx
          !zdif = slope1 * (dist_rx**2+(dist_rx*LeftOrientation(2))**2)**0.5
          !zdif = zdif - (Slope2 * dist_rx*LeftOrientation(2) )
          zdif = slope1 * dist_rx*m1!(dist_rx**2+(dist_rx*m)**2)**0.5
          zdif = zdif + (Slope2 * dist_rx *m2)


          ! point at which everything is flat on second trunk
          slopedist = RotFrontLR(2)-left_yorigin-flatdistance
          extra_slope = zdif/slopedist

          IF(ydif < slopedist) THEN
            extra_z = ydif * extra_slope
          ELSE
            extra_z = zdif
          END IF

          elevation = rz0 - Slope2 * ydif - extra_z
!#endif
          !elevation = -900.0_dp
        ELSE

          elevation = (z0 - slope1*y)
          extra_z=0.0_dp
          zdif=0.0_dp

        END IF

      ELSE ! right

        inside = PointInPolygon(RightPoly, (/x,y/), 0.1_dp, .TRUE.)

        IF(inside) THEN

          elevation = -900.0_dp

!#if 0
          !rotatated coord system with origin as start of second trunk
          NodeHolder(1) = Mesh % Nodes % x(i)
          NodeHolder(2) = Mesh % Nodes % y(i)
          NodeHolder(3) = 0.0_dp

          NodeHolder = MATMUL(RRotationMatrix,NodeHolder)

          rx = rightnodewidth - (NodeHolder(2) - RotFrontRL(1))
          ry = NodeHolder(3) - right_yorigin

          ! distance along divide between trunks
          !addition = (rx**2+(rx*RightOrientation(2))**2)**0.5
          !now need to put thisd 'additon' distance along the divide vector
          !m = (cross(2)-right_lengthtobif)/(cross(1)-width)
          !addition = -1*m*rx


          ! distance in flow direction since divide
          !ydif = ry - rx*LeftOrientation(2)
          dist_rx = rightnodewidth-rx
          m1 = (cross(2)-right_lengthtobif) / rightnodewidth
          m2 = (right_right_min-right_left_min) / rightnodewidth
          ydif =  ry - dist_rx * m2 !- left_y
          addition = m1*(rx)


          ! height at divide
          rz0 = (z0 - slope1 * (right_lengthtobif+addition))


          ! additional slope twisting. Zdiff is the additional heightloss required to flatten out
          ! second trunk. It is distance above lowest point of ice divide minus height loss from extra length
          zdif = slope1 * dist_rx*m1!(dist_rx**2+(dist_rx*m)**2)**0.5
          zdif = zdif + (Slope2 * dist_rx *m2)

          ! point at which everything is flat on second trunk
          slopedist = RotFrontRR(2)-right_yorigin-flatdistance
          extra_slope = zdif/slopedist

          IF(ydif < slopedist) THEN
            extra_z = ydif * extra_slope
          ELSE
            extra_z = zdif
          END IF

          elevation = rz0 - Slope2 * ydif - extra_z
!#endif

        ELSE

          elevation = (z0 - slope1*y)
          extra_z=0.0_dp
          zdif=0.0_dp
        END IF

      END IF
      
    END IF

    Values(Perm(i)) = elevation
  END DO

END SUBROUTINE CreateBed

SUBROUTINE CreateSurface( Model, Solver, dt, TransientSimulation )
  USE CalvingGeometry
  USE MainUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!-------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: NullSolver => NULL()
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  INTEGER :: i,NNodes,dummyint,ierr
  INTEGER, POINTER :: Perm(:), RightPerm(:)=>NULL()
  REAL(KIND=dp), POINTER :: Values(:)
  REAL(KIND=dp) :: FrontOrientation(3),RotationMatrix(3,3),UnRotationMatrix(3,3),Slope1,Slope2,&
    elevation,x,y,z0,width,rx,ry,lengthtobif,NodeHolder(3),rz0,FrontRight(2),FrontLeft(2),RFrontLeft(2),&
    RFrontRight(2),addition,nodewidth,yorigin,ydif,extra_slope,extra_z,slopedist,zdif,flatdistance,dist_rx
  LOGICAL :: Found
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SolverName = "CreateSurface"

  Mesh => Model % Mesh
  NNodes = Mesh % NumberOfNodes

  Values => Solver % Variable % Values
  Perm => Solver % Variable % Perm

  Params => Solver % Values

  ! read in variables
  Slope1 = ListGetConstReal(Params, "Main trunk slope",Found, UnfoundFatal=.TRUE.)
  Slope2 = ListGetConstReal(Params, "Second trunk slope",Found, UnfoundFatal=.TRUE.)
  Width = ListGetConstReal(Params, "Main trunk width",Found, UnfoundFatal=.TRUE.)
  z0 = ListGetConstReal(Params, "z0",Found, UnfoundFatal=.TRUE.)
  flatdistance = ListGetConstReal(Params, "Second trunk flat distance",Found, UnfoundFatal=.TRUE.)

  ! get front orientation
  FrontOrientation = GetFrontOrientation(Model)
  RotationMatrix = ComputeRotationMatrix(FrontOrientation)
  UnRotationMatrix = TRANSPOSE(RotationMatrix)

  ! get front corners
  CALL GetFrontCorners(Model, Solver, FrontLeft, FrontRight)

  ! orientate front corners
  NodeHolder(1:2) = FrontRight
  NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(RotationMatrix,NodeHolder)
  RFrontRight = NodeHolder(2:3)

  NodeHolder(1:2) = FrontLeft
  NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(RotationMatrix,NodeHolder)
  RFrontLeft = NodeHolder(2:3)
  nodewidth = ABS(RFrontLeft(1)-RFrontRight(1))

  !get min point of z origin in rotated coords
  CALL MakePermUsingMask( Model, Solver, Mesh, "Right Sidewall Mask", &
    .FALSE., RightPerm, dummyint)
  yorigin = HUGE(1.0_dp); lengthtobif = HUGE(1.0_dp)
  DO i=1,NNodes
    IF(RightPerm(i) == 0) CYCLE
    x = Mesh % Nodes % x(i)
    y = Mesh % Nodes % y(i)

    NodeHolder(1) = Mesh % Nodes % x(i)
    NodeHolder(2) = Mesh % Nodes % y(i)
    NodeHolder(3) = 0.0_dp

    NodeHolder = MATMUL(RotationMatrix,NodeHolder)
    ry = NodeHolder(3)

    if(ry < yorigin) yorigin = ry
    IF(y < lengthtobif) lengthtobif = y
  END DO

  ! parallel comm
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,yorigin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,lengthtobif, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)

  ! loop through nodes assigning values
  DO i=1, NNodes
    x = Mesh % Nodes % x(i)
    y = Mesh % Nodes % y(i)

    IF(x <= Width) THEN ! on main trunk
      elevation = (z0 - slope1*y)
      extra_z=0.0_dp
      zdif=0.0_dp
    ELSE ! on right trunk

      !rotatated coord system with origin as start of second trunk
      NodeHolder(1) = Mesh % Nodes % x(i)
      NodeHolder(2) = Mesh % Nodes % y(i)
      NodeHolder(3) = 0.0_dp

      NodeHolder = MATMUL(RotationMatrix,NodeHolder)

      rx = nodewidth - (NodeHolder(2) - RFrontLeft(1))
      ry = NodeHolder(3) - yorigin

      ! distance along divide between trunks
      addition = (rx**2+(rx*FrontOrientation(2))**2)**0.5

      ! distance in flow direction since divide
      ydif = ry - rx*FrontOrientation(2)

      ! height at divide
      rz0 = (z0 - slope1 * (lengthtobif+addition))

      ! additional slope twisting. Zdiff is the additional heightloss required to flatten out
      ! second trunk. It is distance above lowest point of ice divide minus height loss from extra length
      dist_rx = nodewidth-rx
      zdif = slope1 * (dist_rx**2+(dist_rx*FrontOrientation(2))**2)**0.5
      zdif = zdif - (Slope2 * dist_rx*FrontOrientation(2) )

      ! point at which everything is flat on second trunk
      slopedist = RFrontLeft(2)-yorigin-flatdistance
      extra_slope = zdif/slopedist

      IF(ydif < slopedist) THEN
        extra_z = ydif * extra_slope
      ELSE
        extra_z = zdif
      END IF

      elevation = rz0 - Slope2 * ydif - extra_z
      
    END IF

    Values(Perm(i)) = elevation

  END DO

END SUBROUTINE CreateSurface

SUBROUTINE CreateFriction( Model, Solver, dt, TransientSimulation )
  USE CalvingGeometry
  USE MainUtils
  USE BedRoutines
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!-------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: NullSolver => NULL()
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  INTEGER :: i,NNodes,dummyint,ierr
  INTEGER, POINTER :: Perm(:), RightPerm(:)=>NULL(),LeftPerm(:)=>NULL()
  REAL(KIND=dp), POINTER :: Values(:)
  REAL(KIND=dp) :: LeftOrientation(3),LRotationMatrix(3,3),LUnRotationMatrix(3,3),&
    RightOrientation(3),RRotationMatrix(3,3),RUnRotationMatrix(3,3),Slope1,Slope2,&
    elevation,x,y,z0,width,rx,ry,NodeHolder(3),rz0,FrontLL(2),FrontLR(2),&
    FrontRL(2),FrontRR(2),RotFrontLL(2),RotFrontLR(2),RotFrontRL(2),RotFrontRR(2),rightnodewidth,&
    addition,leftnodewidth,ydif,extra_slope,extra_z,slopedist,zdif,flatdistance,dist_rx,&
    left_lengthtobif,right_lengthtobif,left_yorigin,right_yorigin,cross(2),LeftPoly(2,5),RightPoly(2,5),&
    m,left_y,LFriction,RFriction,MFriction,left_right_min,left_left_min,right_left_min,right_right_min
  LOGICAL :: Found,inside
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,LeftLeftName,LeftRightName,LeftFrontName,RightLeftName,&
    RightRightName,RightFrontName

  SolverName = "CreateBed"

  Mesh => Model % Mesh
  NNodes = Mesh % NumberOfNodes

  Values => Solver % Variable % Values
  Perm => Solver % Variable % Perm

  Params => Solver % Values

  Width = ListGetConstReal(Params, "Main trunk width",Found, UnfoundFatal=.TRUE.)

  ! read in variables
  LeftLeftName = ListGetString(Params,"Left Channel Left Name", UnfoundFatal=.TRUE.)
  LeftRightName = ListGetString(Params,"Left Channel Right Name", UnfoundFatal=.TRUE.)
  LeftFrontName = ListGetString(Params,"Left Channel Front Name", UnfoundFatal=.TRUE.)

  RightLeftName = ListGetString(Params,"Right Channel Left Name", UnfoundFatal=.TRUE.)
  RightRightName = ListGetString(Params,"Right Channel Right Name", UnfoundFatal=.TRUE.)
  RightFrontName = ListGetString(Params,"Right Channel Front Name", UnfoundFatal=.TRUE.)

  LFriction = ListGetConstReal(Params, "Left channel friction",Found, UnfoundFatal=.TRUE.)
  RFriction = ListGetConstReal(Params, "Right channel friction",Found, UnfoundFatal=.TRUE.)
  MFriction = ListGetConstReal(Params, "Main channel friction",Found, UnfoundFatal=.TRUE.)

  ! get left orientation
  LeftOrientation = GetFrontOrientationCopy(Model, LeftFrontName, LeftLeftName, LeftRightName)
  LRotationMatrix = ComputeRotationMatrix(LeftOrientation)
  LUnRotationMatrix = TRANSPOSE(LRotationMatrix)

  ! get right orientation
  RightOrientation = GetFrontOrientationCopy(Model, RightFrontName, RightLeftName, RightRightName)
  RRotationMatrix = ComputeRotationMatrix(RightOrientation)
  RUnRotationMatrix = TRANSPOSE(RRotationMatrix)

  ! get front corners
  CALL GetFrontCornersCopy(Model, Solver, LeftFrontName, LeftLeftName, LeftRightName, FrontLL, FrontLR)
  CALL GetFrontCornersCopy(Model, Solver, RightFrontName, RightLeftName, RightRightName, FrontRL, FrontRR)

  ! orientate front corners
  NodeHolder(1:2) = FrontLR; NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(LRotationMatrix,NodeHolder)
  RotFrontLR = NodeHolder(2:3)

  NodeHolder(1:2) = FrontLL; NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(LRotationMatrix,NodeHolder)
  RotFrontLL = NodeHolder(2:3)
  leftnodewidth = ABS(RotFrontLL(1)-RotFrontLR(1))

  NodeHolder(1:2) = FrontRR; NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(RRotationMatrix,NodeHolder)
  RotFrontRR = NodeHolder(2:3)

  NodeHolder(1:2) = FrontRL; NodeHolder(3) = 0.0_dp
  NodeHolder = MATMUL(RRotationMatrix,NodeHolder)
  RotFrontRL = NodeHolder(2:3)
  rightnodewidth = ABS(RotFrontRL(1)-RotFrontRR(1))

  !get min point of z origin in rotated coords

  !left first
  CALL MakePermUsingMask( Model, Solver, Mesh, LeftRightName, &
    .FALSE., RightPerm, dummyint)
  CALL MakePermUsingMask( Model, Solver, Mesh, LeftLeftName, &
    .FALSE., LeftPerm, dummyint)

    left_lengthtobif = HUGE(1.0_dp);
  left_left_min = HUGE(1.0_dp); left_right_min = HUGE(1.0_dp)
  DO i=1,NNodes
    IF(RightPerm(i) == 0 .AND. LeftPerm(i) == 0) CYCLE
    x = Mesh % Nodes % x(i)
    y = Mesh % Nodes % y(i)

    NodeHolder(1) = Mesh % Nodes % x(i)
    NodeHolder(2) = Mesh % Nodes % y(i)
    NodeHolder(3) = 0.0_dp

    NodeHolder = MATMUL(LRotationMatrix,NodeHolder)
    ry = NodeHolder(3)

    IF(LeftPerm(i) > 0) THEN
      IF(ry < left_left_min) left_left_min = ry
    END IF

    IF(RightPerm(i) > 0) THEN
      IF(ry < left_right_min) left_right_min = ry
    END IF

    IF(LeftPerm(i) == 0) CYCLE
    IF(y < left_lengthtobif) left_lengthtobif = y
  END DO

  ! parallel comm
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,left_lengthtobif, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,left_right_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,left_left_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)

  left_yorigin = left_right_min!MIN(left_left_min, left_right_min)

  DEALLOCATE(LeftPerm,RightPerm)

  !now right
  CALL MakePermUsingMask( Model, Solver, Mesh, RightRightName, &
    .FALSE., RightPerm, dummyint)
  CALL MakePermUsingMask( Model, Solver, Mesh, RightLeftName, &
    .FALSE., LeftPerm, dummyint)

  right_lengthtobif = HUGE(1.0_dp); cross = HUGE(1.0_dp)
  right_left_min = HUGE(1.0_dp); right_right_min = HUGE(1.0_dp)
  DO i=1,NNodes
    IF(RightPerm(i) == 0 .AND. LeftPerm(i) == 0) CYCLE
    x = Mesh % Nodes % x(i)
    y = Mesh % Nodes % y(i)

    NodeHolder(1) = Mesh % Nodes % x(i)
    NodeHolder(2) = Mesh % Nodes % y(i)
    NodeHolder(3) = 0.0_dp

    NodeHolder = MATMUL(RRotationMatrix,NodeHolder)
    ry = NodeHolder(3)

    IF(LeftPerm(i) > 0) THEN
      IF(ry < right_left_min) right_left_min = ry
    END IF

    IF(RightPerm(i) > 0) THEN
      IF(ry < right_right_min) right_right_min = ry
    END IF

    IF(x<cross(1) .AND. LeftPerm(i) > 0) THEN
      cross(1) = x
      cross(2) = y
    END IF

    IF(RightPerm(i) == 0) CYCLE
    IF(y < right_lengthtobif) right_lengthtobif = y
  END DO

  ! parallel comm
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,right_lengthtobif, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,cross, 2, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,right_right_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,right_left_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)

  right_yorigin = right_left_min!MIN(right_left_min, right_right_min)

  DEALLOCATE(LeftPerm,RightPerm)

  IF(ABS(left_lengthtobif-right_lengthtobif) > AEPS) CALL FATAL(SolverName, 'Bifurcation start at different trunk lengths')

  ! create polygon of left birfucation
  LeftPoly(1,1) = 0.0_dp; LeftPoly(2,1) = left_lengthtobif
  LeftPoly(:,2) = FrontLL
  LeftPoly(:,3) = FrontLR
  LeftPoly(:,4) = cross
  LeftPoly(:,5) = LeftPoly(:,1)

  ! create polygon of left birfucation
  RightPoly(1,1) = width; RightPoly(2,1) = right_lengthtobif
  RightPoly(:,2) = FrontRR
  RightPoly(:,3) = FrontRL
  RightPoly(:,4) = cross
  RightPoly(:,5) = RightPoly(:,1)

  ! loop through nodes assigning values
  DO i=1, NNodes
    x = Mesh % Nodes % x(i)
    y = Mesh % Nodes % y(i)

    inside = PointInPolygon(LeftPoly, (/x,y/), 0.1_dp, .TRUE.)

    IF(inside) THEN
      Values(Perm(i)) = LFriction
    ELSE 
      inside = PointInPolygon(RightPoly, (/x,y/), 0.1_dp, .TRUE.)

      IF(inside) THEN
        Values(Perm(i)) = RFriction
      ELSE
        Values(Perm(i)) = MFriction
      END IF
    END IF
  
  END DO

END SUBROUTINE CreateFriction

SUBROUTINE GetBifurcation( Model, Solver, dt, TransientSimulation )
  USE CalvingGeometry
  USE MainUtils
  USE BedRoutines
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!-------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: VeloVar
  INTEGER, POINTER :: TopPerm(:)=>NULL(),SidePerm(:)=>NULL(),VeloPerm(:),NodeIndexes(:)
  REAL(kind=dp), POINTER :: VeloValues(:)
  INTEGER :: i,j,n,NBulk,NBdry,VeloDOFs,dummyint,BCTag,Idx1,Idx2,counter,nloc,ierr
  REAL(kind=dp) :: Velo1,Velo2,dist,location,sumvelo,x1,x2,locations(10)
  LOGICAL :: Found,ThisBC,OnSide(10),FileCreated=.FALSE.
  INTEGER, ALLOCATABLE :: disps(:),P_nloc(:)
  REAL(kind=dp), ALLOCATABLE :: P_locations(:)
  LOGICAL, ALLOCATABLE :: P_onside(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: SideMaskName,TopMaskName,VeloVarName,SolverName,FileName

  SAVE :: FileCreated

  SolverName = "GetBifurcation"

  Params => Solver % Values
  Mesh => Model % Mesh
  NBulk = Mesh % NumberOfBulkElements
  NBdry = Mesh % NumberOfBoundaryElements

  BCTag = ListGetInteger(Params,"Internal Boundary Number", UnfoundFatal=.TRUE.)
  SideMaskName = ListGetString(Params,"Main Channel Name", UnfoundFatal=.TRUE.)

  !Get the flow solution
  VeloVarName = ListGetString(Params, "Flow Solution Variable Name", Found)
  IF(.NOT. Found) THEN
    CALL Info(SolverName, "Flow Solution Variable Name not found, assuming 'Flow Solution'")
    VeloVarName = "Flow Solution"
  END IF
  VeloVar => VariableGet(Mesh % Variables, VeloVarName, .TRUE., UnfoundFatal=.TRUE.)
  VeloValues => VeloVar % Values
  VeloPerm => VeloVar % Perm
  VeloDOFs = VeloVar % DOFs

  TopMaskName = "Top Surface Mask"
  CALL MakePermUsingMask( Model, Solver, Mesh, TopMaskName, &
      .FALSE., TopPerm, dummyint)

  CALL MakePermUsingMask( Model, Solver, Mesh, SideMaskName, &
      .FALSE., SidePerm, dummyint)

  nloc = 0; OnSide = .FALSE.
  !loop through elements finding one that has both +ve and -ve x velocity
  DO i=NBulk+1, NBulk+NBdry
    IF(Mesh % Elements(i) % BoundaryInfo % Constraint /= BCTag) CYCLE ! not internal bc

    NodeIndexes => Mesh % Elements(i) % NodeIndexes

    IF(COUNT(TopPerm(NodeIndexes) > 0) /= 2) CYCLE ! not two surface nodes present

    n = Mesh % Elements(i) % TYPE % NumberOfNodes

    counter=0
    DO j=1,n
      IF(TopPerm(NodeIndexes(j)) /= 0) CYCLE ! not top node
      counter = counter + 1
      !NodeVelo(1) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 1)
      IF(counter == 1) Idx1 = NodeIndexes(j)
      IF(counter == 2) Idx2 = NodeIndexes(j)
    END DO

    Velo1 = VeloValues(((VeloPerm(Idx1)-1)*VeloDOFs) + 1)
    Velo2 = VeloValues(((VeloPerm(Idx2)-1)*VeloDOFs) + 1)

    IF(Velo1 > 0 .AND. Velo2 > 0) CYCLE
    IF(Velo1 < 0 .AND. Velo2 < 0) CYCLE

    x1 = Mesh % Nodes % x(Idx1)
    x2 = Mesh % Nodes % x(Idx2)

    IF(Velo1 == 0 .AND. SidePerm(Idx1) == 0) THEN
      location = x1

      IF(Velo2 == 0) CALL FATAL(SolverName, 'Element has multiple zero velocities')

    ELSE IF(Velo2 == 0 .AND. SidePerm(Idx2) == 0) THEN
      location = x2

    ELSE

      sumvelo = ABS(Velo1) + ABS(Velo2)
      dist = ABS(x1 - x2)

      IF(x1 > x2) THEN
        location = x1 - ( ABS(Velo1) / sumvelo ) * dist
      ELSE
        location = x1 + ( ABS(Velo1) / sumvelo ) * dist
      END IF

    END IF

    nloc = nloc + 1

    IF(nloc > 10) CALL FATAL(SolverName, 'Too many potential bifurcation points')

    locations(nloc) = location
    
    IF(SidePerm(Idx1) > 0 .OR. SidePerm(Idx2) > 0) OnSide(nloc) = .TRUE.
  END DO

  ALLOCATE(P_nloc(ParEnv % PEs), disps(ParEnv % PEs))

  CALL MPI_AllGather(nloc, 1, MPI_INTEGER, P_nloc, &
            1, MPI_INTEGER, ELMER_COMM_WORLD, ierr)

  disps(1) = 0
  DO i=2, ParEnv % PEs
    disps(i) = disps(i-1) + P_nloc(i-1)
  END DO

  n = SUM(P_nloc)
  ALLOCATE(P_locations(SUM(P_nloc)), P_onside(SUM(P_nloc)))

  !gather corner nodenums
  CALL MPI_AllGatherV(locations(:nloc), nloc, MPI_DOUBLE_PRECISION, P_locations, P_nloc, &
            disps, MPI_DOUBLE_PRECISION, ELMER_COMM_WORLD, ierr)

    !gather corner nodenums
  CALL MPI_AllGatherV(onside(:nloc), nloc, MPI_LOGICAL, P_onside, P_nloc, &
            disps, MPI_LOGICAL, ELMER_COMM_WORLD, ierr)

  IF(ParEnv % MyPE == 0) THEN
    Filename = ListGetString(Params,"Output File Name", Found)
    IF(.NOT. Found) THEN
      CALL WARN(SolverName, 'Output file name not given so using Flux.txt')
      Filename = "Bifurcation.txt"
    END IF

    ! write to file
    IF(FileCreated) THEN
      OPEN( 56, FILE=filename, STATUS='UNKNOWN', POSITION='APPEND')
    ELSE
      OPEN( 56, FILE=filename, STATUS='UNKNOWN')
      WRITE(56, *) "Bifurcation output file for BC ID", BcTag
      WRITE(56, '(A)') "TimeStep, Time, x-coordinate"
    END IF

    !Write out the left and rightmost points
    IF(COUNT(P_ONside) > 0) THEN
      DO i=1,SUM(P_nloc)
        IF(P_Onside(i)) CYCLE
        WRITE(56, *) GetTimestep(), GetTime(), P_locations(i)
      END DO
    ELSE
      DO i=1,SUM(P_nloc)
        WRITE(56, *) GetTimestep(), GetTime(), P_locations(i), 'side!'
      END DO
    END IF

    CLOSE(56)
  END IF

  FileCreated = .TRUE.

END SUBROUTINE GetBifurcation

SUBROUTINE GlacierVolume( Model, Solver, dt, TransientSimulation )

  USE MainUtils
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!-------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: ElevVar
  TYPE(Element_t), POINTER :: Element
  INTEGER, POINTER :: ElevPerm(:),NodeIndexes(:),TopPerm(:)=>NULL()
  REAL(KIND=dp), POINTER :: ElevValues(:)
  INTEGER :: i,NBulk,NBdry,BCId,n,ierr,NNodes,TopCount
  REAL(KIND=dp) :: TotalVolume,Volume,ElemArea,ElemElev,depths,MeanDepth,Zs(10),TotalArea
  LOGICAL :: Found,ThisBC,FileCreated
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,FileName

  SolverName = "GlacierVolume"

  Mesh => Model % Mesh
  Params => Solver % Values
  Nbulk = Mesh % NumberOFBulkElements
  NBdry = Mesh % NumberOFBoundaryElements
  NNodes = MEsh % NumberOfNOdes

  ElevVar => VariableGet( Mesh % Variables, "Elevation", .TRUE. )
  IF(.NOT. ASSOCIATED(ElevVar)) THEN
    CALL Fatal(SolverName, "Couldn't find elevation variable")
  ELSE
    ElevValues => ElevVar % Values
    ElevPerm => ElevVar % Perm
  END IF
  
  DO i=1,Model % NumberOfBCs
    ThisBC = ListGetLogical(Model % BCs(i) % Values,"Top Surface Mask",Found)
    IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
    BCId =  Model % BCs(i) % Tag
    EXIT
  END DO

  CALL MakePermUsingMask( Model, Solver, Mesh, "Top Surface Mask", &
      .FALSE., TopPerm, TopCount)

  depths=0.0_dp
  DO i=i, NNodes
    IF(TopPerm(i) == 0) CYCLE
    depths = depths + ElevValues(ElevPerm(i))
  END DO

  MeanDepth = depths/TopCount

  IF(.NOT. Found) CALL FATAL(SolverName, "Unable to find 'Top Surface Mask'")

  TotalVolume = 0.0_dp; TotalArea = 0.0_dp
  DO i=NBulk+1, NBulk+NBdry
    Element => Mesh % Elements(i)
    IF(Element % BoundaryInfo % Constraint /= BCId) CYCLE ! not target BC

    n = Element % TYPE % NumberOfNodes
    NodeIndexes => Element % NodeIndexes

    Zs(1:n) = Mesh % Nodes % z(NodeIndexes)
    Mesh % Nodes % z(NodeIndexes) = 0.0_dp 

    ElemArea = ElementArea( Solver % Mesh, Element, n )

    Mesh % Nodes % z(NodeIndexes) = Zs(1:n)

    ElemElev = SUM(ElevValues(ElevPerm(NodeIndexes))) / n

    Volume = ElemElev*ElemArea
    
    TotalArea = TotalArea + ElemArea
    TotalVolume = TotalVolume + Volume
  END DO



  CALL MPI_ALLREDUCE(MPI_IN_PLACE,TotalVolume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,TotalArea, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr)

  MeanDepth = TotalVolume/TotalArea

  IF(ParEnv % MYPE == 0) THEN
    Filename = ListGetString(Params,"Output File Name", Found)
    IF(.NOT. Found) THEN
      CALL WARN(SolverName, 'Output file name not given so using Flux.txt')
      Filename = "Volume.txt"
    END IF

    INQUIRE( file=trim(filename), exist=FileCreated )

    ! write to file
    IF(FileCreated .AND. GetTimestep() > 1) THEN
      OPEN( 11, FILE=filename, STATUS='UNKNOWN', POSITION='APPEND')
    ELSE
      OPEN( 11, FILE=filename, STATUS='UNKNOWN')
      WRITE(11, *) "Total volume output for entire domain"
      WRITE(11, '(A)') "TimeStep, Time, Volume, MeanDepth, Area"
    END IF

    !Write out the left and rightmost points
    WRITE(11, *) GetTimestep(), GetTime(), TotalVolume, MeanDepth, TotalArea

    CLOSE(11)
  END IF

END SUBROUTINE GlacierVolume

RECURSIVE SUBROUTINE DummySolver( Model,Solver,Timestep,TransientSimulation )
  USE DefUtils
 
  IMPLICIT NONE


  !------------------------------------------------------------------------------
  !    External variables
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), Pointer :: BC
  TYPE(Variable_t), POINTER :: Var
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: VarPerm(:)
  INTEGER :: VarDOFs, i, j, k, N, t
  REAL(KIND=dp), POINTER :: VarValues(:)
  LOGICAL :: GotIt

  PRINT *,"DummySolver"
  PRINT *,"***********************************"

  Var => Solver % Variable
  IF (ASSOCIATED(Var)) THEN
     VarPerm => Var % Perm
     VarDOFs =  Var % DOFs
     VarValues => Var % Values
  ELSE
     CALL FATAL('DummySolver','No Variable associated')
  END IF
  k=0
 ! DO i = 1,Model % NumberOfNodes
 !    DO j=1,VarDOFs
 !       k = k + 1
 !       VarValues(VarDOFs*(VarPerm(i) - 1)+j) = k
 !    END DO
 ! END DO
 ! VarValues = 0.0_dp
  DO t=1, Solver % Mesh % NumberOfBoundaryElements
      ! get element information
      Element => GetBoundaryElement(t)
      IF ( .NOT.ActiveBoundaryElement() ) CYCLE
      BC => GetBC()
      n = GetElementNOFNodes()
      VarValues(VarPerm(Element % NodeIndexes)) = ListGetReal(BC,TRIM(Solver % Variable % Name),n,Element % NodeIndexes,GotIt)
   END DO
  
END SUBROUTINE DummySolver