module types
   implicit none

   ! FOR THE GNU FORTRAN COMPILER
   integer,parameter,public :: l_1 = 2  ! standart Logical type
   integer,parameter,public :: i_2 = 2  ! 16 bits integer
   integer,parameter,public :: i_4 = 4  ! 32 bits integer
   integer,parameter,public :: r_4 = 4  ! 32 bits float
   integer,parameter,public :: r_8 = 8  ! 64 bits float

end module types

program allocation !module to test allocation logic module of LPJmL-Fire
    
    use types
    use iso_fortran_env, only : output_unit

    !VARIABLES [INPUT] - Determinadas arbitrariamente

    integer(i_4),parameter :: npls = 20
    integer(i_4),parameter :: nseg = 20 ! segment number to bisection loop
    real(r_8),parameter :: pi   =  3.14159265
    real(r_8),parameter :: xacc =  0.1     !x-axis precision threshold for the allocation solution
    real(r_8),parameter :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution
    integer(i_4),parameter :: ntl=365
    integer, parameter :: stdout = output_unit


    integer(i_4) :: pls
    integer(i_4) :: xmin, xmax
    real(r_8) :: allom1 = 100
    real(r_8) :: allom2 = 40.0
    real(r_8) :: allom3 = 0.5
    real(r_8) :: latosa = 8.e3
    real(r_8) :: reinickerp = 1.6
    real(r_8) :: ltor = 0.77302587552347657 !leaf:root from Philip
    real(r_8), dimension(npls) :: dwood   !in g/cm3 but must be in g/m2
    real(r_8), dimension(npls) :: sla !m2/g
    real(r_8), dimension(npls) :: nind !m2
    real(r_8), dimension(npls) :: bminc !kgC/m2 - !total biomass increment this year for area
    real(r_8), dimension(npls) :: height
    

    real(r_8), dimension(npls) :: lm_ind !in kgC/m2 but use in gC/ind [transformation below]
    real(r_8), dimension(npls) :: sm_ind !in kgC/m2 but use in gC/ind [transformation below]
    real(r_8), dimension(npls) :: hm_ind !in kgC/m2 but use in gC/ind [transformation below]
    real(r_8), dimension(npls) :: rm_ind !in kgC/m2 but use in gC/ind [transformation below]
    real(r_8), dimension(npls) :: cw_ind !to calculate sap and heartwood (in kgC/m2 but use in gC/ind)
    real(r_8), dimension(npls) :: litter_ag_fast
    real(r_8), dimension(npls) :: litter_ag_slow
    real(r_8), dimension(npls) :: litter_bg

    !Local Variables
    real(r_8),dimension(npls) :: bminc_ind !gC/ind - individual total biomass increment this year 
    real(r_8),dimension(npls) :: lm2rm          !ratio of leafmass to fine rootmass
    real(r_8),dimension(npls) :: lminc_ind      !individual leafmass increment this year
    real(r_8),dimension(npls) :: rminc_ind      !individual fineroot mass increment this year
    real(r_8),dimension(npls) :: lminc_ind_min  !min leafmass increment to maintain current sapwood
    real(r_8),dimension(npls) :: rminc_ind_min  !min rootmass increment to support new leafmass
    real(r_8),dimension(npls) :: sap_xsa        !cross sectional area of sapwood  
    real(r_8),dimension(npls) :: sminc_ind      !individual sapmass increment this year

    real(r_8),dimension(npls) :: lm !leaf mass
    real(r_8),dimension(npls) :: sm !sapwood mass
    real(r_8),dimension(npls) :: hm !heartwood mass
    real(r_8),dimension(npls) :: rm !root mass

    real(r_8),dimension(npls) :: x1             !working vars in bisection
    real(r_8),dimension(npls) :: x2
    real(r_8),dimension(npls) :: rtbis
    real(r_8),dimension(npls) :: dx
    real(r_8),dimension(npls) :: xmid
    real(r_8),dimension(npls) :: root1, root2, root3
    real(r_8),dimension(npls) :: sign
    real(r_8) :: wooddens = 2.e5
    logical  :: normal

    real(r_8),dimension(npls) :: fx1
    real(r_8),dimension(npls) :: fmid

    real(r_8),dimension(npls) :: lm1     !allometric leafmass requirement (leafmass req'd to keep sapwood alive; gC ind-1)

    integer :: i


    !Arrays with values to some variables (generic values)
    dwood =(/0.74,0.73,0.59,0.52,0.41,0.44,0.86,0.42,0.64,0.69,0.92,&
    &0.60,0.36,0.99,0.59,0.52,0.41,0.44,0.86,0.42/) !atenção para a unidade
    bminc=(/2.15,2.,2.18,2.6,2.5,1.8,2.3,2.,1.8,2.84,2.25,3.,2.2,1.7,&
    &1.18,2.6,3.5,2.8,3.3,2./)
    sla=(/0.002,0.018,0.009,0.023,0.013,0.039,0.040,0.0028,0.0025,&
    &0.027,0.032,0.007,0.013,0.025,0.002,0.008,0.004,0.016,0.023,0.015/)
    nind=(/1.,2.,8.,6.,5.,9.,3.,4.,7.,1.,2.,8.,5.,3.,6.,4.,5.,8.,9.,3./)
    height=(/5.,9.,15.,10.9,11.5,18.9,12.6,2.5,14.9,22.5,28.7,23.6,&
    &28.8,19.6,13.3,27.6,29.5,21.6,30.,2./)
    lm_ind=(/2.15,2.,1.18,1.6,1.5,1.8,0.3,2.,0.8,.84,0.25,1.,0.2,1.7,&
    &1.18,1.6,1.5,1.8,0.3,2./)
    cw_ind=(/7.,12.,7.2,8.3,8.8,9.7,7.5,11.5,10.,8.6,7.3,10.3,6.8,9.9,&
    &5.3,9.2,15.,12.6,10.7,11.4/)
    rm_ind=(/0.63,0.8,0.9,0.5,1.3,0.9,0.4,1.0,0.56,0.87,0.33,0.97,0.31,&
    &0.55,0.2,0.8,0.4,0.66,0.23,1.5/)

    !-----------------------------------------------------------------


    do pls = 1,npls

        lm(pls) = (lm_ind(pls)/nind(pls))*1.D3 !PRECISA COLOCAR VALORES INICIAIS
        sm(pls) = ((cw_ind(pls)*0.05)/nind(pls))*1.D3
        hm(pls) = ((cw_ind(pls)*0.95)/nind(pls))*1.D3
        rm(pls) = (rm_ind(pls)/nind(pls))*1.D3

        ! print*, 'LM=', lm(pls),'SM',sm(pls),'HM', hm(pls),'RM', rm(pls)
        
        bminc_ind(pls) = (bminc(pls)/nind(pls))*1.D3

        ! ====== TREE ALLOCATION ======

        ! lm1(pls) = latosa*sm(pls)/(dwood*height(pls)*sla(pls))  !allometric leaf mass requirement *****ATENÇÃO*****

        lm1(pls) = 1000.0  !valor arbitrario colocado para rever a unidade do dwood

        lminc_ind_min(pls) = lm1(pls) - lm(pls)  !eqn (27)
    
        !calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
        !May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.

        rminc_ind_min(pls) = lm1(pls) / ltor - rm(pls)      !eqn (30)


        if (rminc_ind_min(pls) .gt. 0. .and. lminc_ind_min(pls) .gt. 0. .and. &
            &(rminc_ind_min(pls) + lminc_ind_min(pls)) .le. bminc_ind(pls)) then

            !Normal allocation (positive increment to all living C compartments)

            normal = .true.

            !Calculation of leaf mass increment (lminc_ind) that satisfies Eqn (22)
            !Since this is normal allocation, we set the lower bound for the leafmass allocation (x1)
            !to its allometric minimum, because it should be able to be fulfilled, i.e.:

            !Start to find root procedure (relate to bisection method)

            x1(pls) = lminc_ind_min(pls)
            x2(pls) = (bminc_ind(pls) - (lm(pls) / ltor - rm(pls))) / (1. + 1. / ltor)
            
            dx(pls) = x2(pls) - x1(pls)

            if (dx(pls) < 0.01) then

                !there seems to be rare cases where lminc_ind_min (x1) is almost equal to x2. In this case,
                !assume that the leafmass increment is equal to the midpoint between the values and skip 
                !the root finding procedure

                lminc_ind(pls) = x1(pls) + 0.5 * dx(pls)

            else
                !Find a root for non-negative lminc_ind, rminc_ind and sminc_ind using Bisection Method (Press et al 1986, p 346)
                !There should be exactly one solution (no proof presented, but Steve has managed one).
                    
                dx(pls) = dx(pls)/nseg

                !! ===== FIND ROOT FUNCTION ===== [**must be a function**]

                pi4 = pi/4
                a1 = 2./allom3
                a2 = 1. + a1
                a3 = allom2**a1


                root1(pls) = a3*((sm(pls)+bminc_ind(pls)-x1(pls)-((lm(pls)+x1(pls))/ltor)+&
                        &rm(pls)+hm(pls))/wooddens)/pi4-((sm(pls)+bminc_ind(pls)-x1(pls)-&
                        &((lm(pls)+x1(pls))/ltor)+rm(pls))/((lm(pls)+x1(pls))*sla(pls)*&
                        &wooddens/latosa))**a2

                ! ======================================================

                !evaluate f(x1) = LHS of eqn (22) at x1

                fx1(pls) = root1(pls)

                !Find approximate location of leftmost root on the interval (x1,x2).
                !Subdivide (x1,x2) into nseg equal segments seeking change in sign of f(xmid) relative to f(x1).

                fmid(pls) = fx1(pls)
                xmid(pls) = x1(pls)

                i = 1

                do

                    xmid(pls) = xmid(pls) + dx(pls)

                    root2(pls) = a3*((sm(pls)+bminc_ind(pls)-xmid(pls)-((lm(pls)+xmid(pls))/ltor)+&
                    &rm(pls)+hm(pls))/wooddens)/pi4-((sm(pls)+bminc_ind(pls)-xmid(pls)-&
                    &((lm(pls)+xmid(pls))/ltor)+rm(pls))/((lm(pls)+xmid(pls))*sla(pls)*&
                    &wooddens/latosa))**a2

                    fmid(pls) = root2(pls)

                    if ((fmid(pls)*fx1(pls)) .le. 0. .or. xmid(pls) .ge. x2(pls)) exit  !sign has changed or we are over the upper bound

                    if (i > 20) write(stdout,*)'first alloc loop flag',i,pls,fmid(pls)*fx1(pls),&
                         &xmid(pls),x1(pls),x2(pls),dx(pls),bminc_ind(pls)

                    if (i > 50) stop 'Too many iterations allocmod'

                    i = i + 1

                end do

                !the interval that brackets zero in f(x) becomes the new bounds for the root search

                x1(pls) = xmid(pls) - dx(pls)
                x2(pls) = xmid(pls)

                !Apply bisection method to find root on the new interval (x1,x2)

                fx1(pls) = root1(pls)

                if (fx1(pls) .ge. 0.) then
                    sign(pls) = -1.
                else
                    sign(pls) =  1.
                end if

                rtbis(pls) = x1(pls)
                dx(pls)    = x2(pls) - x1(pls)

                !Bisection loop: search iterates on value of xmid until xmid lies within xacc of the root,
                !i.e. until |xmid-x| < xacc where f(x) = 0. the final value of xmid with be the leafmass increment

                i = 1

                do 

                    dx(pls)   = 0.5 * dx(pls)
                    xmid(pls) = rtbis(pls) + dx(pls)

                    !calculate fmid = f(xmid) [eqn (22)]

                    root3(pls) = a3*((sm(pls)+bminc_ind(pls)-xmid(pls)-((lm(pls)+xmid(pls))/ltor)+&
                    &rm(pls)+hm(pls))/wooddens)/pi4-((sm(pls)+bminc_ind(pls)-xmid(pls)-&
                    &((lm(pls)+xmid(pls))/ltor)+rm(pls))/((lm(pls)+xmid(pls))*sla(pls)*&
                    &wooddens/latosa))**a2

                    fmid(pls) = root3(pls)

                    if (fmid(pls) * sign(pls) .le. 0.) rtbis(pls) = xmid(pls)

                    if (dx(pls) < xacc .or. abs(fmid(pls)) <= yacc) exit

                    if (i > 20) write(stdout,*)'second alloc loop flag',i,pls,dx(pls),abs(fmid(pls))
                    if (i > 50) stop 'Too many iterations allocmod'

                    i = i + 1

                end do

                !Now rtbis contains numerical solution for lminc_ind given eqn (22)

                lminc_ind(pls) = rtbis(pls)

            endif

            !Calculate increments in other compartments using allometry relationships

            rminc_ind(pls) = (lm(pls) + lminc_ind(pls)) / ltor - rm(pls)       !eqn (9)

            sminc_ind(pls) = bminc_ind(pls) - rminc_ind(pls) - lminc_ind(pls)  !eqn (1)

            ! print*, 'LEAF_INC (gC/ind)', (lminc_ind(pls)/1.D3), 'ROOT_INC (gC/ind)', (rminc_ind(pls)/1.D3),&
            ! & 'SAP_INC(gC/ind)', (sminc_ind(pls)/1.D3), pls

        else 

            !Abnormal allocation: reduction in some C compartment(s) to satisfy allometry
            
            normal = .false.

            !Attempt to distribute this year's production among leaves and roots only

            lminc_ind(pls) = (bminc_ind(pls)-lm(pls)/ltor+rm(pls))/(1.+1./ltor)  !eqn (33)


            if (lminc_ind(pls) > 0.) then

                !Positive allocation to leafmass

                rminc_ind(pls) = bminc_ind(pls) - lminc_ind(pls)  !eqn (31)
                
                !Add killed roots (if any) to below-ground litter

                if (rminc_ind(pls) < 0.) then

                    lminc_ind(pls) = bminc_ind(pls)
                    rminc_ind(pls) = (lm(pls) + lminc_ind(pls)) / ltor - rm(pls)

                    litter_bg(pls) = litter_bg(pls) + abs(rminc_ind(pls)) * nind(pls)

                end if
                
                i = 1

            else

                !Negative allocation to leaf mass

                rminc_ind(pls) = bminc_ind(pls)
                lminc_ind(pls) = (rm(pls) + rminc_ind(pls)) * ltor - lm(pls)  !from eqn (9)

                !Add killed leaves to litter

                litter_ag_fast(pls) = litter_ag_fast(pls) + abs(lminc_ind(pls)) * nind(pls)
                
                i = 2

            endif

            !Calculate sminc_ind (must be negative)
      
            sminc_ind(pls) = (lm(pls) + lminc_ind(pls)) * sla(pls) /&
            & latosa * 2.e5 * height(pls) - sm(pls)  !eqn (35)

            !Convert killed sapwood to heartwood

            hm(pls) = hm(pls) + abs(sminc_ind(pls))


        endif !normal/abnormal allocation

        !Increment C compartments - OUTPUT FINAL (kgC/m²)

        lm_ind(pls) = ((lm(pls) + lminc_ind(pls))*nind(pls))/1.D3
        rm_ind(pls) = ((rm(pls) + rminc_ind(pls))*nind(pls))/1.D3 !PQ TA IGUAL AS FOLHAS AFF
        sm_ind(pls) = ((sm(pls) + sminc_ind(pls))*nind(pls))/1.D3
        hm_ind(pls) = (hm(pls)*nind(pls))/1.D3

        print*, 'LM', lm_ind(pls), 'RM', rm_ind(pls), 'SM', sm_ind(pls), 'HM', hm_ind(pls), pls

    enddo

    


end program allocation