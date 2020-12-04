module allometry
	use ISO_FORTRAN_ENV, only: REAL32, REAL64, REAL128
	use constants
	use allocation

	implicit none

	contains

	function diam (dw, pi,delta_sapwood) result (diameter) !!!!!precisamos rever todas as unidades
		real(REAL64) :: diameter 
		real(REAL64) :: dw 
		real(REAL64) :: pi
		real(REAL64) :: delta_sapwood !Deve ser o valor de carbono total no sapwood.

		diameter = ((4 + (delta_sapwood)) / ((dw) * pi * 40)) **&
		&(1 / (2 + 0.5))
		
	end function diam

	function area_crown (k_allom1, krp, diameter)  result (crown_area)
		real(REAL64) :: diameter
		real(REAL64) :: k_allom1
		real(REAL64) :: krp
		real(REAL64) :: crown_area

		crown_area = k_allom1*(diameter**krp)

	end function area_crown

	function tree_height (k_allom2, k_allom3, diameter) result (height)
		real(REAL64) :: diameter
		real(REAL64) :: k_allom2
		real(REAL64) :: k_allom3
		real(REAL64) :: height

		height = k_allom2*(diameter**k_allom3)

	end function tree_height


end module allometry