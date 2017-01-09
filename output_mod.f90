MODULE output

USE MPI_module
USE variables
USE set_pressure

IMPLICIT NONE

INTEGER(MPI_OFFSET_KIND), DIMENSION(:, :), ALLOCATABLE :: displ
INTEGER(MPI_OFFSET_KIND), DIMENSION(:),    ALLOCATABLE :: loc_displ
INTEGER, DIMENSION(:, :, :), ALLOCATABLE               :: shapes
INTEGER, DIMENSION(:, :),    ALLOCATABLE               :: loc_shapes
INTEGER, DIMENSION(:), ALLOCATABLE                     :: MPI_block_type


CONTAINS


  SUBROUTINE set_output
    ! The only purpose of this subroutine is to define the displacements necessary
    ! to each process to understand where it can begin writing its own stuff.

    IMPLICIT NONE

    INTEGER                            :: i, j, ierr
    INTEGER                            :: dblsz, ss
    INTEGER(MPI_ADDRESS_KIND)          :: i1, i2
    INTEGER, DIMENSION(ndims)          :: array_of_sizes
    INTEGER, DIMENSION(ndims)          :: array_of_subsizes
    INTEGER, DIMENSION(ndims)          :: array_of_starts

    ! displacements
    ALLOCATE(displ(n_flow_variables, 0:nprocs-1))
    ALLOCATE(loc_displ(n_flow_variables))
    ! shapes
    ALLOCATE(shapes(n_flow_variables, ndims, 0:nprocs-1))
    ALLOCATE(loc_shapes(n_flow_variables, ndims))
    ! exchange type
    ALLOCATE(MPI_block_type(n_flow_variables))

    ! debug-purpose initialization
    displ = -99
    shapes = -99

    ! size of a double precision variable
    dblsz = SIZEOF(uvwp(1)%values(uvwp(1)%b(1, 1), &
                                  uvwp(1)%b(2, 1), &
                                  uvwp(1)%b(3, 1)))

    !!!!!!!!!! Sizes in bytes of each block !!!!!!!!!!
    ! velocity components
    DO i = 1, n_flow_variables-1
      loc_displ(i) = dblsz*(uvwp(i)%b_bc(1, 2)-uvwp(i)%b_bc(1, 1)+1)*&
                           (uvwp(i)%b_bc(2, 2)-uvwp(i)%b_bc(2, 1)+1)*&
                           (uvwp(i)%b_bc(3, 2)-uvwp(i)%b_bc(3, 1)+1)
      loc_shapes(i, :) = [uvwp(i)%b_bc(1, 2)-uvwp(i)%b_bc(1, 1)+1, &
                          uvwp(i)%b_bc(2, 2)-uvwp(i)%b_bc(2, 1)+1, &
                          uvwp(i)%b_bc(3, 2)-uvwp(i)%b_bc(3, 1)+1]
    END DO
    ! pressure
    i = n_flow_variables
    !loc_displ(i) = SIZEOF(uvwp(i)%values)
    loc_displ(i) = SIZEOF(p)
    loc_shapes(i, :) = SHAPE(p)

    ! every process needs to know this info (not really, but that's how it works now)
    CALL MPI_ALLGATHER(loc_displ, n_flow_variables, MPI_LONG, &
                       displ,     n_flow_variables, MPI_LONG, &
                       MPI_COMM_WORLD, ierr)
    ! from chunk sizes to displacements of file pointer
    for_each_process: DO i = nprocs-1, 1, -1
      for_each_variable: DO j = 1, n_flow_variables
        displ(j, i) = SUM(displ(j, 0:i-1))
      END DO for_each_variable
    END DO for_each_process
    i = 0
    displ(:, i) = 0


    ! Giving all shape data to master process
    CALL MPI_GATHER(loc_shapes, n_flow_variables*ndims, MPI_INTEGER, &
                    shapes,     n_flow_variables*ndims, MPI_INTEGER, &
                    0, MPI_COMM_WORLD, ierr)


    !!!!!!!!!! type definition !!!!!!!!!!
    for_each_vel_comp: DO i = 1, n_flow_variables-1

      array_of_sizes = SHAPE(uvwp(i)%values)
      array_of_starts = [0, 0, 0]
      array_of_subsizes = [uvwp(i)%b_bc(1, 2)-uvwp(i)%b_bc(1, 1)+1, &
                           uvwp(i)%b_bc(2, 2)-uvwp(i)%b_bc(2, 1)+1, &
                           uvwp(i)%b_bc(3, 2)-uvwp(i)%b_bc(3, 1)+1]
      CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
             & array_of_sizes, array_of_subsizes, array_of_starts, &
             & MPI_ORDER_FORTRAN, &
             & MPI_DOUBLE_PRECISION, MPI_block_type(i), &
             & ierr)
      CALL MPI_TYPE_COMMIT(MPI_block_type(i), ierr)

    END DO for_each_vel_comp

    i = n_flow_variables
    !array_of_sizes = SHAPE(uvwp(i)%values)
    array_of_sizes = SHAPE(p)
    array_of_starts = [0, 0, 0]
    !array_of_subsizes = SHAPE(uvwp(i)%values)
    array_of_subsizes = SHAPE(p)
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, &
           & array_of_sizes, array_of_subsizes, array_of_starts, &
           & MPI_ORDER_FORTRAN, &
           & MPI_DOUBLE_PRECISION, MPI_block_type(i), &
           & ierr)
    CALL MPI_TYPE_COMMIT(MPI_block_type(i), ierr)

  END SUBROUTINE set_output



  SUBROUTINE raw_out
    ! This subroutine writes a file of raw data by means of MPI-I/O tools. Being
    ! written in binary, it requires additional info to be of any use.

    IMPLICIT NONE

    INTEGER :: ierr, i
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER :: file

    !TODO: dare in output un unico file, accedendo così una sola volta all'hard disk

    !!!!!!!!!! u !!!!!!!!!!
    i = 1
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, 'u', MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file, ierr)
    CALL MPI_FILE_SEEK(file, displ(i, myid), MPI_SEEK_SET, ierr)
    CALL MPI_FILE_WRITE(file, uvwp(i)%values(uvwp(i)%b_bc(1, 1),  &
                                             uvwp(i)%b_bc(2, 1),  &
                                             uvwp(i)%b_bc(3, 1)), &
                                             1, MPI_block_type(i), status, ierr)
    CALL MPI_FILE_CLOSE(file, ierr)

    !!!!!!!!!! v !!!!!!!!!!
    i = 2
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, 'v', MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file, ierr)
    CALL MPI_FILE_SEEK(file, displ(i, myid), MPI_SEEK_SET, ierr)
    CALL MPI_FILE_WRITE(file, uvwp(i)%values(uvwp(i)%b_bc(1, 1),  &
                                             uvwp(i)%b_bc(2, 1),  &
                                             uvwp(i)%b_bc(3, 1)), &
                                             1, MPI_block_type(i), status, ierr)
    CALL MPI_FILE_CLOSE(file, ierr)

    !!!!!!!!!! w !!!!!!!!!!
    i = 3
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, 'w', MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file, ierr)
    CALL MPI_FILE_SEEK(file, displ(i, myid), MPI_SEEK_SET, ierr)
    CALL MPI_FILE_WRITE(file, uvwp(i)%values(uvwp(i)%b_bc(1, 1),  &
                                             uvwp(i)%b_bc(2, 1),  &
                                             uvwp(i)%b_bc(3, 1)), &
                                             1, MPI_block_type(i), status, ierr)
    CALL MPI_FILE_CLOSE(file, ierr)

    !!!!!!!!!! p !!!!!!!!!!
    i = 4
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, 'p', MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file, ierr)
    CALL MPI_FILE_SEEK(file, displ(i, myid), MPI_SEEK_SET, ierr)
    !CALL MPI_FILE_WRITE(file, uvwp(i)%values, 1, MPI_block_type(i), status, ierr)
    CALL MPI_FILE_WRITE(file, p, 1, MPI_block_type(i), status, ierr)
    CALL MPI_FILE_CLOSE(file, ierr)


    CALL output_read_aid


  END SUBROUTINE raw_out



  SUBROUTINE output_read_aid

    USE grids, ONLY: L

    IMPLICIT NONE

    INTEGER                :: iostat
    INTEGER, PARAMETER     :: dev = 0
    INTEGER                :: i

    IF (myid==0) THEN

      OPEN(UNIT=dev, FILE='outputaid.txt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=iostat)

      WRITE(dev, *) 'Cartesian topology'
      WRITE(dev, *) dims
      WRITE(dev, *) ''

      WRITE(dev, *) 'Domain dimensions'
      WRITE(dev, *) L
      WRITE(dev, *) ''

      WRITE(dev, *) 'u shapes'
      DO i = 0, nprocs-1
        WRITE(dev, *) shapes(1, :, i)
      END DO
      WRITE(dev, *) ''

      WRITE(dev, *) 'v shapes'
      DO i = 0, nprocs-1
        WRITE(dev, *) shapes(2, :, i)
      END DO
      WRITE(dev, *) ''

      WRITE(dev, *) 'w shapes'
      DO i = 0, nprocs-1
        WRITE(dev, *) shapes(3, :, i)
      END DO
      WRITE(dev, *) ''

      WRITE(dev, *) 'p shapes'
      DO i = 0, nprocs-1
        WRITE(dev, *) shapes(4, :, i)
      END DO
      WRITE(dev, *) ''

    END IF

  END SUBROUTINE output_read_aid




END MODULE output
