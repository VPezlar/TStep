MODULE Arnoldi

USE accuracy
USE setup


    implicit none

contains

    subroutine arnoldi_iter(A, n, m, v, H)
        ! Implements the Arnoldi Eigenvalue Algorithm
        ! A: Input matrix (n x n)
        ! n: Dimension of the matrix
        ! m: Krylov subspace size
        ! v: Output basis vectors (n x (m+1))
        ! H: Output Hessenberg matrix (m+1 x m)

        REAL(rk), dimension(n, n),   intent(in)  :: A
        INTEGER(ik),                 intent(in)  :: n, m
        REAL(rk), dimension(n, m+1), intent(out) :: v
        REAL(rk), dimension(m+1, m), intent(out) :: H

        REAL(rk), dimension(n) :: w
        REAL(rk)               :: norm_w
        INTEGER(ik)                     :: j, i

        ! Initialize H to zeros
        H = 0.0_rk

        ! Step 1: Normalize initial random vector v1 (assuming v(:,1) is already set and non-zero)
        ! For testing, we'll assume v(:,1) is provided and needs normalization.
        ! If you want a random v1, you'd generate it here.

        norm_w = sqrt(sum(v(:,1)**2))
        if (norm_w > tiny(1.0_rk)) then
            v(:,1) = v(:,1) / norm_w
        else
            ! Handle case where initial vector is zero
            write(*,*) "Error: Initial vector v1 is zero. Cannot normalize."
            return
        end if


        ! Step 2: for j = 1 to m do
        do j = 1, m
            ! Step 3: w = Av_j
            w = matmul(A, v(:,j))

            ! Step 4: for i = 1 to j do
            do i = 1, j
                ! Step 5: h_i,j = v_i^T w
                H(i,j) = dot_product(v(:,i), w)

                ! Step 6: w = w - h_i,j v_i (Orthogonalization)
                w = w - H(i,j) * v(:,i)
            end do ! end for i

            ! Step 8: h_j+1,j = ||w||
            norm_w = sqrt(sum(w**2))
            H(j+1,j) = norm_w

            ! Step 9: if h_j+1,j != 0 then
            if (H(j+1,j) > tiny(1.0_rk)) then
                ! Step 10: v_j+1 = w/h_j+1,j
                v(:,j+1) = w / H(j+1,j)
            else
                ! If norm_w is zero, breakdown occurred.
                ! We can stop early or handle this case as needed.
                write(*,*) "Arnoldi breakdown at j =", j, ". ||w|| is zero."
                ! For now, we'll break the loop.
                exit
            end if
        end do ! end for j

    end subroutine arnoldi_iter

end module Arnoldi