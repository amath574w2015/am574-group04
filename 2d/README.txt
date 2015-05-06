
This is very incomplete....

rpn2euq3.f, rpt2euq.f, setaux.f were all copied from old 4.x code.  I fixed
the calling sequences, but have not re-ordered all the indices.

bc2quad.f also from old 4.x code, implemented solid wall boundaries for a
mapped grid and inflow boundary condition for nozzle.  Rather than modifying
this it might be easiest to incorporate these changes into a version of the
new library routine $CLAW/classic/src/2d/bc2.f90.  

For an axisymmetric nozzle, you will also need to add a source term.

