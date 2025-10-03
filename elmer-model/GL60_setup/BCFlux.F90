SUBROUTINE BCFlux( Model, Solver, dt, TransientSimulation )
  USE CalvingGeometry
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
  TYPE(Element_t), POINTER :: Element
  TYPE(Variable_t), POINTER :: VeloVar
  TYPE(Nodes_t) :: Nodes
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: i,j,NBulk,NBdry,BCId,n,idx,ierr
  REAL(KIND=dp) :: Area,Flux,ElemArea,Normal(3),NodeVelo(3),ElemFlux
  LOGICAL :: Found,FirstTime,FileCreated=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,VeloVarName,FileName

  !SAVE :: FileCreated

  SolverName = "BCFlux"

  Mesh => Model % Mesh
  NBulk = Mesh % NumberOfBulkElements
  NBdry = Mesh % NumberOfBoundaryElements
  Params => Solver % Values

  ! get id tag
  BCId = GetInteger(Params, 'Target Boundary', Found )
  IF(.NOT. Found) CALL FATAL(SolverName, "Unable to find 'Target Boundary'")

  !Get the flow solution
  VeloVarName = ListGetString(Params, "Flow Solution Variable Name", Found)
  IF(.NOT. Found) THEN
    CALL Info(SolverName, "Flow Solution Variable Name not found, assuming 'Flow Solution'")
    VeloVarName = "Flow Solution"
  END IF
  VeloVar => VariableGet(Mesh % Variables, VeloVarName, .TRUE., UnfoundFatal=.TRUE.)

  Area = 0.0_dp
  Flux = 0.0_dp

  !really if applied properly should only enter on one bc but lets double check
  n=0
  DO i=NBulk+1, NBulk+NBdry
    Element => Mesh % Elements(i)
    IF(Element % BoundaryInfo % Constraint /= BCId) CYCLE ! not target BC

    FirstTime = (n /= Element % TYPE % NumberOfNodes) ! different number of nodes?

    n = Element % TYPE % NumberOfNodes
    ElemArea = ElementArea( Solver % Mesh, Element, n )

    NodeIndexes => Element % NodeIndexes

    IF(FirstTime) THEN
      IF(ASSOCIATED(Nodes % x)) DEALLOCATE(Nodes % x)
      IF(ASSOCIATED(Nodes % y)) DEALLOCATE(Nodes % y)
      IF(ASSOCIATED(Nodes % z)) DEALLOCATE(Nodes % z)
      ALLOCATE(Nodes % x(n), Nodes % y(n), Nodes % z(n))
    END IF

    Nodes % x = Mesh % Nodes % x(NodeIndexes)
    Nodes % y = Mesh % Nodes % y(NodeIndexes)
    Nodes % z = Mesh % Nodes % z(NodeIndexes)

    Normal = NormalVector(Element, Nodes)

    NodeVelo = 0.0_dp
    DO j=1,n
      idx = NodeIndexes(j)
      NodeVelo(1) = NodeVelo(1) + VeloVar % Values(((VeloVar % Perm(idx)-1)*VeloVar % DOFs) + 1)
      NodeVelo(2) = NodeVelo(2) + VeloVar % Values(((VeloVar % Perm(idx)-1)*VeloVar % DOFs) + 2)
      NodeVelo(3) = NodeVelo(3) + VeloVar % Values(((VeloVar % Perm(idx)-1)*VeloVar % DOFs) + 3)
    END DO
    NodeVelo = NodeVelo/n

    ElemFlux = DOT_PRODUCT(NodeVelo,Normal)
    ElemFlux = ElemFlux * ElemArea

    Flux = Flux + ElemFlux
    Area = Area + ElemArea
  END DO

  CALL MPI_ALLREDUCE(MPI_IN_PLACE,Flux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,Area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr)

  IF(ParEnv % MYPE == 0) THEN
    Filename = ListGetString(Params,"Output File Name", Found)
    IF(.NOT. Found) THEN
      CALL WARN(SolverName, 'Output file name not given so using Flux.txt')
      Filename = "Flux.txt"
    END IF

    INQUIRE( file=trim(filename), exist=FileCreated )

    ! write to file
    IF(FileCreated .AND. GetTimestep() > 1) THEN
      OPEN( 64, FILE=filename, STATUS='UNKNOWN', POSITION='APPEND')
    ELSE
      OPEN( 64, FILE=filename, STATUS='UNKNOWN')
      WRITE(64, *) "Flux output file for BC ID", BcId
      WRITE(64, '(A)') "TimeStep, Time, Area, Total Flux, Mean Flux"
    END IF

    !Write out the left and rightmost points
    WRITE(64, *) GetTimestep(), GetTime(), Area, Flux, Flux/Area

    CLOSE(64)
  END IF
  
  IF(ASSOCIATED(Nodes % x)) DEALLOCATE(Nodes % x)
  IF(ASSOCIATED(Nodes % y)) DEALLOCATE(Nodes % y)
  IF(ASSOCIATED(Nodes % z)) DEALLOCATE(Nodes % z)

  !FileCreated = .TRUE.

END SUBROUTINE BCFlux