MODULE solve_pressure
! The present module contains the following procedures
!    SUBROUTINE Solve_p
!    SUBROUTINE Comm_p
!    SUBROUTINE Res_check (tol,norm, done, it_num)
!    SUBROUTINE Out_diag
!    SUBROUTINE P_iter
!    SUBROUTINE Rel_vec(omega)

    USE grids
    USE MPI_module
    USE set_pressure
    USE mpi

    IMPLICIT NONE

    INTEGER :: it_num = 0, M
    REAL :: norm
    LOGICAL :: done = .FALSE.

    REAL, ALLOCATABLE, DIMENSION(:) :: omega

CONTAINS

    !!!!!!!!!!Routine chiamata nel MAIN
    SUBROUTINE Solve_p

        IMPLICIT NONE

        IF (It_type .eq. 'SRJ') CALL Rel_vec(omega)

        ! the iteration cycle must be entered at least once
        done = .FALSE.
        it_num = 0
        DO WHILE (done .EQV. .FALSE.)

            !Comunicazione delle p di faccia.
            CALL Comm_p

            !Residuo PRIMA della nuova iterazione
            CALL Res_check (tol, norm, done, it_num)

            !Nuova p
            CALL P_iter

        END DO

    END SUBROUTINE Solve_p


    SUBROUTINE Comm_p

        IMPLICIT NONE

        INTEGER :: ierr,i
        INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat

        DO i = 1,ndims

            CALL MPI_SENDRECV(p_sen(i)%p_p, SIZE(p_sen(i)%p_p), MPI_DOUBLE_PRECISION, idp(i), 11, &
                            & p_rec(i)%p_m, SIZE(p_rec(i)%p_m), MPI_DOUBLE_PRECISION, idm(i), 11, procs_grid, stat, ierr)

            CALL MPI_SENDRECV(p_sen(i)%p_m, SIZE(p_sen(i)%p_m), MPI_DOUBLE_PRECISION, idm(i), 22, &
                            & p_rec(i)%p_p, SIZE(p_rec(i)%p_p), MPI_DOUBLE_PRECISION, idp(i), 22, procs_grid, stat, ierr)

        END DO

    END SUBROUTINE Comm_p


    SUBROUTINE Res_check (tol,norm, done, it_num)

        IMPLICIT NONE

        REAL, INTENT(IN) :: tol
        LOGICAL, INTENT(OUT) :: done
        INTEGER, INTENT (INOUT) :: it_num
        REAL, INTENT(OUT) :: norm
        REAL :: norm_loc
        INTEGER :: ierr

        it_num = it_num + 1

        CALL Out_diag!(res)

        !Definisci b, OD
        res = OD+weigC_grid*p-b

        norm_loc = MAXVAL(ABS(res))
        CALL MPI_ALLREDUCE(norm_loc, norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, procs_grid, ierr)
        IF (norm<= tol) done = .TRUE.

        !IF (myid == 0) WRITE(*,*) norm, it_num

    END SUBROUTINE Res_check


    SUBROUTINE Out_diag
    ! TODO: cercare di ridurre tutte queste operazioni utilizzando la routine di
    ! moltiplicazione CDS x vettore

        IMPLICIT NONE

        !real, ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) :: OD_dummy
        INTEGER, DIMENSION(3) :: Up, Lo
        INTEGER :: i, j, k

        DO k = indv(3,1), indv(3,2)

            DO j = indv(2,1), indv(2,2)

                DO i = indv(1,1), indv(1,2)

                    Up(1:1) = MIN(Lapl(1)%ud, ovl_rec(1,2)-i)
                    Lo(1:1) = MIN(Lapl(1)%ld, i-ovl_rec(1,1))
                    Up(2:2) = MIN(Lapl(2)%ud, ovl_rec(2,2)-j)
                    Lo(2:2) = MIN(Lapl(2)%ld, j-ovl_rec(2,1))
                    Up(3:3) = MIN(Lapl(3)%ud, ovl_rec(3,2)-k)
                    Lo(3:3) = MIN(Lapl(3)%ld, k-ovl_rec(3,1))

                    OD(i,j,k)  =  SUM(Lapl(1)%matrix(i, 1:Up(1))*uvwp(4)%values(i+1:i+Up(1),j,k)) &
                    & + SUM(Lapl(1)%matrix(i,-Lo(1):-1)*uvwp(4)%values(i-Lo(1):i-1,j,k)) &
                    !
                    & + SUM(Lapl(2)%matrix(j, 1:Up(2))*uvwp(4)%values(i,j+1:j+Up(2),k)) &
                    & + SUM(Lapl(2)%matrix(j,-Lo(2):-1)*uvwp(4)%values(i, j-Lo(2):j-1,k)) &
                    !
                    & + SUM(Lapl(3)%matrix(k, 1:Up(3))*uvwp(4)%values(i,j, k+1:k+Up(3)))  &
                    & + SUM(Lapl(3)%matrix(k,-Lo(3):-1)*uvwp(4)%values(i, j, k-Lo(3):k-1))

                END DO

            END DO

        END DO

    END SUBROUTINE Out_diag


    SUBROUTINE P_iter

        IMPLICIT NONE

        INTEGER :: j

        iteration_method: SELECT CASE(It_type)

            CASE('J')

                p = p - res/weigC_grid

            CASE('SRJ')

                j = MOD(it_num, M)
                IF (j == 0) j = M
                p = p - omega(j)*res/weigC_grid

        END SELECT iteration_method

    END SUBROUTINE P_iter


    SUBROUTINE Rel_vec(omega)

        IMPLICIT NONE

        TYPE ind_type
            INTEGER, ALLOCATABLE, DIMENSION(:) :: row
        END TYPE ind_type

        INTEGER :: Lev
        TYPE(ind_type), ALLOCATABLE, DIMENSION(:) :: ind_vec
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: omega
        REAL, ALLOCATABLE, DIMENSION(:) ::om_appo
        INTEGER, ALLOCATABLE, DIMENSION(:) :: q
        INTEGER :: i, Nmax

        Nmax = MAXVAL([Ntot(1), Ntot(2), Ntot(3)])

        SELECT CASE(Nmax)

            CASE (16:31)

                Lev = 5
                M = 43
                ALLOCATE(omega(M),om_appo(Lev))
                om_appo(1) = 88.190
                om_appo(2) = 30.122
                om_appo(3) = 6.8843
                om_appo(4) = 1.6008
                om_appo(5) = 0.58003

                ALLOCATE(q(Lev-1), ind_vec(Lev-1))
                q(1) = 1
                q(2) = 2
                q(3) = 5
                q(4) = 12

                DO i = 1,Lev-1
                ALLOCATE(ind_vec(i)%row(q(i)))
                END DO
                ind_vec(1)%row = 1
                ind_vec(2)%row = [2,24]
                ind_vec(3)%row = [3,11,20,27,35]
                ind_vec(4)%row = [(4+i*3, i = 0,11)]

            CASE (32:63)
                Lev = 5
                M = 76
                ALLOCATE(omega(M),om_appo(Lev))
                om_appo(1) = 330.57
                om_appo(2) = 82.172
                om_appo(3) = 13.441
                om_appo(4) = 2.2402
                om_appo(5) = 0.60810

                ALLOCATE(q(Lev-1), ind_vec(Lev-1))
                q(1) = 1
                q(2) = 2
                q(3) = 7
                q(4) = 20

                DO i = 1,Lev-1
                ALLOCATE(ind_vec(i)%row(q(i)))
                END DO
                ind_vec(1)%row = 1
                ind_vec(2)%row = [2,41]
                ind_vec(3)%row = [3,14,23,33,44,53,63]
                ind_vec(4)%row = [(4+i*3, i = 0,19)]

            CASE (64:127)

                Lev = 15
                M = 301
                ALLOCATE(omega(M),om_appo(Lev))
                om_appo(1) = 1604.55
                om_appo(2) = 1236.6
                om_appo(3) = 777.72
                om_appo(4) = 429.57
                om_appo(5) = 220.699
                om_appo(6) = 109.268
                om_appo(7) = 53.1653
                om_appo(8) = 25.7023
                om_appo(9) = 12.4395
                om_appo(10) = 6.06839
                om_appo(11) = 3.77684
                om_appo(12) = 2.26342
                om_appo(13) = 1.17188
                om_appo(14) = 0.697364
                om_appo(15) = 0.519746


                ALLOCATE(q(Lev-1), ind_vec(Lev-1))
                q(1) = 1
                q(2) = 1
                q(3) = 1
                q(4) = 2
                q(5) = 2
                q(6) = 4
                q(7) = 6
                q(8) = 9
                q(9) = 13
                q(10) = 20
                q(11) = 7
                q(12) = 36
                q(13) = 51
                q(14) = 68


                DO i = 1,Lev-1
                    ALLOCATE(ind_vec(i)%row(q(i)))
                END DO

                ind_vec(1)%row = 1
                ind_vec(2)%row = 67
                ind_vec(3)%row = 100
                ind_vec(4)%row = [138,164]
                ind_vec(5)%row = [47, 182]
                ind_vec(6)%row = [33,126, 195, 209]
                ind_vec(7)%row = [22, 88, 113, 218, 226, 236]
                ind_vec(8)%row = [16, 59, 76, 144, 155, 177, 204, 243, 248]
                ind_vec(9)%row = [9, 41, 54, 95, 109, 122, 134, 173, 252, 257, 260, 265, 268]
                ind_vec(10)%row = [13, 30, 37, 64, 73, 84, 92, 118, 152, 159, 168, 187, 192,201,215,222,232,241,271,273]
                ind_vec(11)%row = [5, 26, 50, 80, 105, 129, 147]
                ind_vec(12)%row = [7, 18, 28, 43, 52, 61, 71, 86, 98, 107, 116, 131, 141, 150, 161, &
                & 171, 180, 184, 190, 198, 208, 212, 220, 229, 234, 239, 250,  &
                & 256, 263, 276, 277, 279, 280, 282, 283, 285]
                ind_vec(13)%row = [3, 12, 15, 21, 25, 34, 36, 39, 46, 48, 57, 58, 66, 69,77, 79, 82,  &
                & 90, 94, 102, 103, 111, 115, 121, 124, 127, 136, 137, 143, 148, &
                & 156, 158, 166, 167, 176, 178, 189, 197, 203, 206, 217, 225, 228, &
                & 245, 247, 255, 262, 270, 286, 287, 288]
                ind_vec(14)%row = [4, 8, 11, 20, 24, 29, 32, 42, 45, 53, 56, 63, 65, 74, 75, 87, 89, 97, &
                & 99, 101, 110, 114, 120, 132, 133, 135, 146, 151, 154, 163, 165, 172, &
                & 175, 183, 186, 188, 194, 196, 200, 202, 211, 214, 216, 221, 224, 231, &
                & 233, 237, 238, 240, 244, 251, 254, 259, 266, 267, 269, 275, 289, 290, &
                & 291, 292, 293, 294, 295, 296, 297, 298]

            CASE(128:255)
                Lev = 15
                M = 593
                ALLOCATE(omega(M),om_appo(Lev))
                om_appo(1) = 6371.830
                om_appo(2) = 4666.190
                om_appo(3) = 2709.050
                om_appo(4) = 1370.560
                om_appo(5) = 646.1340
                om_appo(6) = 294.6220
                om_appo(7) = 132.3670
                om_appo(8) = 59.13710
                om_appo(9) = 26.41590
                om_appo(10) = 11.85290
                om_appo(11) = 5.942440
                om_appo(12) = 2.897170
                om_appo(13) = 1.364490
                om_appo(14) = 0.7437780
                om_appo(15) = 0.5238630


                ALLOCATE(q(Lev-1), ind_vec(Lev-1))
                q(1) = 1
                q(2) = 1
                q(3) = 1
                q(4) = 2
                q(5) = 3
                q(6) = 5
                q(7) = 8
                q(8) = 12
                q(9) = 19
                q(10) = 29
                q(11) = 30
                q(12) = 68
                q(13) = 101
                q(14) = 141

                DO i = 1,Lev-1
                    ALLOCATE(ind_vec(i)%row(q(i)))
                END DO

                ind_vec(1)%row = 1
                ind_vec(2)%row = 172
                ind_vec(3)%row = 239
                ind_vec(4)%row = [115, 307]
                ind_vec(5)%row = [74, 351, 377]
                ind_vec(6)%row = [50, 20, 283, 406, 427]
                ind_vec(7)%row = [32, 135, 157, 265, 442, 457, 467, 479]
                ind_vec(8)%row = [22, 87, 103, 196, 229, 253, 328, 339, 369, 395, 493, 503]
                ind_vec(9)%row = [15, 58, 67, 129, 148, 166, 186, 222, 291, 299, 320, 421, &
                & 451, 512, 516, 523, 528, 533, 540]
                ind_vec(10)%row = [11, 41, 46, 84, 98, 109, 122, 142, 182, 210, 217, 247,  &
                & 260, 271, 277, 317, 347, 357, 364, 385, 390, 403, 417, &
                & 437, 473, 485, 499, 544, 548]
                ind_vec(11)%row = [7, 28, 38, 64, 79, 95, 126, 151, 162, 179, 199, 203, 236, &
                & 245, 287, 302, 313, 335, 343, 381, 410, 433, 447, 463, &
                & 489, 509, 537, 552, 554, 556]
                ind_vec(12)%row = [5, 18, 25, 36, 54, 57, 62, 72, 77, 92, 106, 113, 118, 120, &
                & 137, 140, 146, 160, 171, 176, 191, 193, 215, 221, 227,  &
                & 232, 235, 257, 263, 269, 274, 281, 294, 297, 311, 325, &
                & 330, 333, 355, 361, 368, 373, 375, 393, 398, 401, 415, &
                & 425, 430, 441, 455, 460, 471, 477, 483, 497, 506, 521, &
                & 525, 531, 559, 560, 562, 564, 565, 567, 569, 571]
                ind_vec(13)%row = [12, 14, 21, 24, 31, 34, 44, 45, 49, 52, 69, 70, 82, 86, 89, &
                & 91, 100, 101, 105, 111, 128, 132, 133, 145, 154, 156, 159, &
                & 165, 169, 175, 184, 188, 189, 201, 205, 209, 212, 214, 225, &
                & 242, 244, 249, 251, 254, 256, 268, 279, 285, 290, 293, 305, &
                & 306, 310, 319, 322, 324, 340, 342, 346, 349, 353, 360, 366, &
                & 383, 386, 388, 392, 408, 412, 414, 420, 423, 435, 439, 445, &
                & 449, 453, 459, 466, 470, 476, 487, 492, 495, 502, 505, 514, &
                & 517, 519, 535, 539, 543, 547, 550, 572, 574, 575, 576, 577, 578]
                ind_vec(14)%row = [9, 10, 17, 20, 29, 30, 39, 40, 43, 56, 60, 61, 63, 65, 66, 75, 76,    &
                & 80, 81, 94, 97, 99, 110, 116, 117, 119, 123, 124, 125, 131, 139, 141, &
                & 144, 149, 150, 153, 164, 168, 174, 178, 181, 183, 195, 197, 198, 200, 208,&
                & 219, 220, 224, 228, 231, 233, 234, 238, 240, 241, 261, 262, 266, 267, 273,&
                & 275, 276, 278, 284, 289, 296, 300, 301, 304, 315, 316, 318, 329, 332, 334,&
                & 337, 338, 345, 356, 359, 363, 365, 371, 372, 376, 379, 380, 382, 396, 397,&
                & 400, 402, 405, 407, 411, 419, 428, 429, 432, 434, 438, 444, 452, 456, 462,&
                & 465, 469, 475, 481, 482, 484, 486, 491, 498, 501, 510, 511, 513, 524, 527,&
                & 529, 530, 534, 538, 546, 555, 558, 563, 570, 579, 580, 581, 582, 583, 584,&
                & 585, 586, 587]

            CASE(256:511)!Manca!!
                Lev = 15
                M = 1154
                ALLOCATE(omega(M),om_appo(Lev))

                om_appo(1) = 25234.4
                om_appo(2) = 17233
                om_appo(3) = 9009.26
                om_appo(4) = 4080.06
                om_appo(5) = 1729.54
                om_appo(6) = 712.786
                om_appo(7) = 290.428
                om_appo(8) = 117.862
                om_appo(9) = 47.8274
                om_appo(10) = 19.4772
                om_appo(11) = 8.32423
                om_appo(12) = 3.5472
                om_appo(13) = 1.54899
                om_appo(14) = 0.785632
                om_appo(15) = 0.527435

                ALLOCATE(q(Lev-1), ind_vec(Lev-1))

                q(1) = 1
                q(2) = 1
                q(3) = 1
                q(4) = 2
                q(5) = 4
                q(6) = 6
                q(7) = 11
                q(8) = 18
                q(9) = 29
                q(10) = 48
                q(11) = 68
                q(12) = 125
                q(13) = 196
                q(14) = 286

                DO i = 1,Lev-1
                    ALLOCATE(ind_vec(i)%row(q(i)))
                END DO

                ind_vec(1)%row = 1
                ind_vec(2)%row = 321
                ind_vec(3)%row = 499
                ind_vec(4)%row = [197, 634]
                ind_vec(5)%row = [123, 442, 709, 752 ]
                ind_vec(6)%row = [75, 275, 394,  598, 815, 844]
                ind_vec(7)%row = [47, 169, 243,  363, 539, 571,  678, 798, 879, 894, 915]
                ind_vec(8)%row = [30, 104, 150,  224, 310, 345,  429, 465, 484,  524, 658, 735,  782, 939, 951, 962, 972, 989]
                ind_vec(9)%row = [19, 65, 94, 140, 187,  214, 266, 290,  299, 386, 409,  420, 564, 591,  613, 625, 650, 700,&
                & 724, 773, 835, 860, 869,  932, 1001,  1006,   1013,   1023,   1030]
                ind_vec(10)%row = [13, 42,   60, 89, 118,    133,    162, 179, 208, 237, 253, 259, 331, 339, 357, 372, 379, 451,&
                & 459,  478, 493, 509, 516, 533, 549, 557, 584, 670, 687, 695, 746, 761, 768, 807, 822, &
                & 829,  900, 908, 926, 982, 1037, 1041, 1047, 1051, 1055, 1060, 1066, 1070]
                ind_vec(11)%row = [8, 26,    38, 56, 80, 85, 110, 114, 146,  157, 176, 194, 205, 230, 234, 282,  286, 295, 307,&
                & 316, 328, 354, 401, 405, 415, 426, 437, 446, 475, 505, 546, 580, 604, 609, 620, 629, 640, &
                & 646, 655, 665, 684, 716, 720, 732, 743, 790, 794, 804, 850, 855, 865, 876, 887, 891, 923, &
                & 948,  957, 968, 978, 998, 1019, 1074, 1077, 1080, 1084, 1086, 1089, 1092]
                ind_vec(12)%row = [6, 17, 23, 35, 51, 54, 69, 72, 78, 98, 101, 107,  128, 131, 138,  143, 155, 167,  174,    &
                & 185, 192, 203, 218, 221, 228, 247, 250, 257, 264, 269, 273, 280, 304, 325, 336, 343, 349, &
                & 352, 367, 370, 376, 383, 389, 392, 400, 424, 435, 456, 463, 470, 472, 487, 491, 496, 503, &
                & 515, 521, 528, 531, 536, 543, 555, 562, 569, 574, 577, 589, 596, 601, 618, 638, 663, 676, &
                & 681, 693, 704, 707, 712, 729, 741, 756, 759, 766, 776, 779, 785, 788, 813, 819, 826, 832, &
                & 840, 842, 848, 874, 885, 906, 913, 918, 921, 936, 943, 945, 965, 986, 993, 996, 1010, &
                & 1016, 1028, 1034, 1045, 1064, 1095, 1097, 1099, 1101, 1103, &
                & 1105, 1107, 1108, 1111, 1112, 1114, 1116]
                ind_vec(13)%row = [4, 11, 15, 22, 29, 33, 40, 44, 46, 50, 62, 64, 68, 83, 88, 91, 93, 97, 113,&
                & 117, 120, 122, 126, 136, 149, 153, 160, 163, 165, 172, 182, 184,  190, 199, 201, 211, 213,    &
                & 216, 226, 236, 240, 242, 245, 256, 263, 278, 288, 291, 293, 297,  301, 303, 312, 314, 318,    &
                & 320, 323, 333, 335, 342, 359, 361, 365, 375, 382, 398, 407, 411,  413, 417, 419, 422, 432,    &
                & 434, 440, 444, 448, 450, 454, 461, 468, 480, 482, 486, 502, 511,  513, 519, 526, 542, 552,    &
                & 553, 560, 567, 583, 587, 594, 607, 611, 615, 616, 623, 627, 631,  633, 637, 644, 648, 652,    &
                & 653, 657, 661, 668, 672, 674, 686, 690, 692, 699, 703, 719, 723,  726, 728, 738, 739, 748,    &
                & 750, 754, 764, 771, 775, 796, 801, 802, 806, 810, 811, 818, 825,  838, 854, 858, 861, 863,    &
                & 867, 871, 872, 882, 884, 893, 897, 899, 902, 904, 911, 929, 930,  935, 942, 954, 956, 960,    &
                & 964, 971, 974, 976, 981, 984, 991, 1004, 1008, 1015, 1022, 1026, 1033,&
                & 1040, 1044, 1053, 1057, 1059, 1062, 1072, 1076, 1082, 1091, &
                & 1118, 1119, 1121, 1122, 1123, 1124, 1125, 1127, 1128, 1129]
                ind_vec(14)%row = [3, 9, 10, 14, 21, 25, 28, 32, 37, 39, 49, 55, 58, 59, 61, 71, 74, 76, 77, 81, &
                & 82, 87, 100,  103, 105, 106,  109, 111, 116, 129, 130, 134, 135, 141, 142, 145, 147, 152, 156,&
                & 159, 170, 171, 177, 178, 181, 189, 193, 196, 198, 206, 207, 210, 220, 223, 225, 231, 232, 233,&
                & 239, 249, 252, 254, 255, 261, 262, 268, 271, 272, 276, 277, 283, 284, 287, 296, 306, 309, 311,&
                & 327, 329, 330, 338, 341, 347, 348, 350, 351, 355, 356, 364, 369, 371, 374, 378, 381, 385, 388,&
                & 390, 393, 396, 397, 403, 404, 410, 416, 427, 428, 431, 439, 443, 453, 458, 460, 466, 467, 474,&
                & 476, 477, 479, 490, 492, 495, 497, 498, 501, 507, 508, 510, 518, 523, 525, 530, 532, 535, 538,&
                & 540, 545, 548, 550, 556, 559, 565, 566, 572, 573, 576, 579, 582, 586, 592, 593, 599, 600, 603,&
                & 606, 610, 621, 622, 628, 636, 642, 643, 647, 660, 664, 667, 671, 679, 680, 683, 689, 697, 698,&
                & 702, 710, 711, 714, 715, 717, 722, 733, 734, 736, 745, 747, 753, 758, 762, 763, 769, 770, 778,&
                & 781, 783, 784, 787, 791, 792, 795, 799, 809, 817, 821, 824, 828, 831, 833, 836, 837, 843, 846,&
                & 847, 851, 853, 857, 866, 877, 878, 881, 889, 890, 896, 909, 910, 916, 917, 920, 922, 925, 927,&
                & 934, 938, 941, 947, 949, 952, 953, 959, 967, 970, 980, 988, 990, 995, 999, 1000, 1002,&
                & 1003, 1012, 1014, 1020, 1021, 1025, 1032, 1036, 1039, 1043,&
                & 1049, 1050, 1052, 1068, 1069, 1071, 1075, 1081, 1087, 1088,&
                & 1094, 1100, 1104, 1110, 1117, 1130, 1131, 1132, 1133, 1134,&
                & 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144]
            CASE(512:1023)
        END SELECT

        omega = om_appo(Lev)
        DO i = 1,Lev-1
            omega(ind_vec(i)%row) = om_appo(i)
        END DO

    END SUBROUTINE Rel_vec

END MODULE solve_pressure
