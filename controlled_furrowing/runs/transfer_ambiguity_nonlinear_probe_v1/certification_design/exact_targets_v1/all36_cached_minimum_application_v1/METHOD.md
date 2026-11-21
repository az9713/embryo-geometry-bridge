# Local minima of 36 fixed cached-quadrature energies

Each candidate state has a separately certified cached-quadrature stationary point, an exact binary preconditioner \(R\), and a whole-cube cached derivative bound. Strict cached-site determinant, hoop, radius, tangent and Jacobian margins keep the scalar normalized cached energy \(E_Q\) analytic on a neighborhood of the declared closed cube. The auxiliary finite-element body geometry evidence serves a separate purpose. It identifies a physically admissible closed displacement cube around the same center.

On each matching cached derivative and geometry cube of radius \(2^{-34}\), let \(H_Q\) be the symmetric Hessian of \(E_Q\) and \(A=RH_Q\). A bound \(\|I-A\|_\infty\le Z_Q<1\) makes \(A\) invertible and positive stable. Every eigenvalue of \(A\) has real part at least \(1-Z_Q\). Thus the square factors \(R,H_Q\) are invertible. With \(S=H_Q^{-1}\) and independently certified \(P=(R+R^{\mathsf T})/2\succ0\), the identity \(AS+SA^{\mathsf T}=2P\) gives

\[
S=\int_0^\infty e^{-tA}\,2P\,e^{-tA^{\mathsf T}}\,dt\succ0 .
\]

Therefore \(H_Q\succ0\). An inverse bound from the Neumann series and the exact row-sum norm of \(R\) gives the quantitative curvature bound \(m=(1-Z_Q)/\|R\|_\infty>0\). Positivity of \(P\) matters. For example, \(R=H_Q=-I\) gives zero derivative defect and a strict local maximum instead.

Use the native root posterior for each state wherever it puts the root strictly inside the radius \(2^{-34}\) cube around the same center, independent of whether its native root cube was smaller, equal or larger. The one state needing an extra localization check is candidate + at small input 02. Its native \(2^{-32}\) posterior is wider than the geometry cube. The independently replayed matching fresh \(2^{-34}\) cached bound and same signed center force give a strict small self-map test. Its stationary point is the same unique root as in the enclosing native cube. This is a scalar use of existing bounds. The independently matched body geometry and cached analytic domain cover the same whole radius \(2^{-34}\) cube.

For an interior stationary point \(u_*\), integrate the Hessian bound along the segment from \(u_*\) to \(u\) in its closed radius \(2^{-34}\) cube. This gives \(E_Q(u)-E_Q(u_*)\ge (m/2)\|u-u_*\|_2^2\), and hence a strict local minimum and the unique minimum on that cube. The two recentered proofs retain the original field and the exact center displacement \(\sigma\). Their frozen observer scores continue to use the original archived \(\xi\) and the shift-inclusive verified error.

Every conclusion applies to a fixed binary cached-quadrature energy and its recorded coordinate scaling. The application asserts no exact-spatial minimum, dynamic loading result, continuum result, experimental compatibility, or new computational mechanics calculation. Promotion of the 36-case application requires separate independent actual scalar readback.
