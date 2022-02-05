# About
This project aims to record the codes and experiences while learning POD(Proper Orthogonal Decomposition).

There are three codes to demonstrate the power and features of the POD, i.e.
- `POD_linear.m` is about a simple linear circuit containing only resistors.
- `POD_RC.m` & `POD_RC2.m` are about a simple linear circuit containing a capacitor and two resistors.
- `POD_transeq.m` is about a transport equation.

The POD is implemented via SVD(Singular Value Decomposition), which can be found in the references for more details.

All codes are written in MATLAB.
# Demo1
The circuit is shown below.

![POD_linear](https://github.com/skywalkerwed/POD_demo/blob/main/POD_linear.svg "Circuit of Demo1")

The equation of this circuit is given by

[1 0 0; -1/R 5/2R -1/R; -1/R -1/R 3/R]\*[u1 u2 u3]'=[u(t) 0 0]'

In the code, R is set to 1 (for simplicity).The corresponding results are given in the following.

![POD_linear_u=e^-t](https://github.com/skywalkerwed/POD_demo/blob/main/Linear_fig1.svg "Result of Demo1")

![POD_linear_u=e^-t_2](https://github.com/skywalkerwed/POD_demo/blob/main/Linear_fig2.svg "Result of Demo1")

![POD_linear_u=sint](https://github.com/skywalkerwed/POD_demo/blob/main/Linear_fig3.svg "Result of Demo1")

Remark:
- All cases have verified the feasibility of POD (if the decomposition is implemented properly).
- The first two cases show that the choice of snapshots seems to have no effect on the results. Well, that's **not true**, and counterexamples will be given later.
- The low-rank approximation transformation matrix is the same in all three cases, i.e. phi = [-0.7741 -0.4764 -0.4168]. This result
can also be obtained by Kirchhoff's current law, that is

  (u1-u2)/R=u2/2R+(u2-u3)/R (for node u2) &&

  (u1-u3)/R+(u2-u3)/R=u3/R (for node u3)

  The above two equation can be rewritten as u1=t && u2=8/13t && u3=7/13t, which has the same meaning of phi.
# Demo2 Part I
The circuit of this demo is shown below. Consider the zero-state case.

![POD_RC](https://github.com/skywalkerwed/POD_demo/blob/main/POD_RC.svg "Circuit of Demo2")

The equation of this circuit is given by

Cu1'+1/R(u1-u2)=0 (u1' is the time derivative of u1) &&

-u1/R+2u2/R-u3/R=0 &&

u3=u(t)

In the code, R is set to 50K and C is set to 10mu. The corresponding results are given in the following.

The first case is done by selecting the proper snapshots.

![POD_RC_u=e^-t](https://github.com/skywalkerwed/POD_demo/blob/main/RC_fig1.svg "Result of Demo2, Part I")

The second and third cases are done by selecting somewhat improper snapshots.

![POD_RC_u=e^-t_2](https://github.com/skywalkerwed/POD_demo/blob/main/RC_fig2.svg "Result of Demo2, Part I")

![POD_RC_u=e^-t_3](https://github.com/skywalkerwed/POD_demo/blob/main/RC_fig3.svg "Result of Demo2, Part I")

Remark:
- The choice of snapshot **does affect** the accuracy of the solution. In the first case, snapshots are selected properly to reflect the R-C circuit.
But in the second and third cases, the system reflected by snapshots is somewhat like a linear system shown in demo1, which is totally far away from the original system.
- Furthermore, it is conceivable that POD is input-dependent, and different inputs **may** result in different reduced models. (Yes, there exist some special cases, such as part II)
# Demo2 Part II
The special case is to reduce the order of the same system, but the low-rank transformation matrix used is obtained in other situations, namly different inputs. 
The corresponding results and, for comparison,  the "correctly" processed results are given below. It's amazing that the accuracy of the two results are so closed.

![POD_RC_u=cost](https://github.com/skywalkerwed/POD_demo/blob/main/RC2_fig1.svg "Result of Demo2, Part II")

![POD_RC_u=cost_special](https://github.com/skywalkerwed/POD_demo/blob/main/RC2_fig2.svg "Result of Demo2, Part II")

Remark:
- The transformation matrix obtained by u=cos(t), called phi, is [0.0874 -0.4102 -0.9078; 0.9087 0.4063 -0.0961]'.
The other obtained by u=e^(-t), called phi2, is [-0.4128 -0.5640 -0.7152; 0.814 0.1234 -0.5673]'. By further elementary column transformation, all transformation matrices can be rewritten as [1 0 -1; 0 1 2].
  This result can also be obtained by Kirchhoff's current law, that's to consider the continuity of the current in the circuit, which has the equation as
  
  (u2-u1)/R=(u3-u2)/R (the current flowing through the first resistor must be equal to the second one)
  
  Let u1 be v, and u3 be u (of course, u3 is connected to the source). The equation above can be can be rewritten as [u1 u2 u3]=u[1 1/2 0]'+v[0 1/2 1]', or [u1 u2 u3]=u[1 0 -1]'+v[0 1 2]'.It means that the phi or phi2 is exactly the basis of the solution space.
- Although the transformation matrices are the "same", the corresponding singular values are quite different, that's the physical interpretation by POD is totally different. Of course, it emphasizes that POD is input-dependent.
# Demo3
In this demo, a simple heat transport is studied, the equation is given as

u'\_t-a^2u''\_xx=0 &&

u(x,0)=u0 &&

u(0,t)=u1 &&

u(l,t)=u2

The result is shown as follow, for better demonstration, it's extended to 2-d with the hot colormap.

![POD_transeq_2](https://github.com/skywalkerwed/POD_demo/blob/main/2-d.gif "Result of Demo3, 2-d")
# References
Pinnau R. (2008) Model Reduction via Proper Orthogonal Decomposition. In: Schilders W.H.A., van der Vorst H.A., Rommes J. (eds) Model Order Reduction: Theory, Research Aspects and Applications. Mathematics in Industry (The European Consortium for Mathematics in Industry), vol 13. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-78841-6_5
# License
The content of this project itself is licensed under the `cc-by-4.0` license, and the codes are licensed under the `mit` license.
