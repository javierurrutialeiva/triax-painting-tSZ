# Triax-Painting-tSZ
Painting algorithm to N-body simulation using a GNFW profile with triaxiality. Equations were partially obtained from [Kim et al 2025](https://arxiv.org/abs/2307.04794) to define elliptical coordinates. y-profile is completelly based on [Battaglia et al 2016](https://iopscience.iop.org/article/10.1088/0004-637X/758/2/75/pdf)
## Projected ellipsoidal coordinates:
Ellipsodial radius is given by:
$$
\zeta^2 = \sqrt{\frac{\hat{x}^2_1}{q_1^2} + \frac{\hat{x}^2_2}{q_2^2} + \hat{x}^3_2}
$$
Prolate if $q_1 = q_2 \leq 1$ and oblate if $q_1 \leq q_2 = 1$. 
#### Euler angles 
$(x^1_{\mathrm{int}}, x^2_{\mathrm{int}}, x^3_{\mathrm{int}})\rightarrow (x^1_{\mathrm{obs}}, x^2_{\mathrm{obs}}, x^3_{\mathrm{obs}})$. There are three angles:

* $\theta$: Angle between $x^3_{\mathrm{int}}$ and $x^3_{\mathrm{obs}}$ aligned with the ellipse mayor-axis.
* $\varphi$ Angle between $x^1_{\mathrm{int}}$ and the line of nodes.
* $\phi$ Angle between $x^1_{\mathrm{obs}}$ and the line of nodes.

The **line of nodes** is the interception of the $x^{1}_{\mathrm{int}}-x^{2}_{\mathrm{int}}$ and $x^{2}_{\mathrm{int}}-x^{1}_{\mathrm{obs}}$ planes. The ellepticity $\epsilon$ of the projected profile becomes $\epsilon = 1 - q_p$, where $q_p$:
$$
q_p = \sqrt{\frac{j + l - \sqrt{(j - l)^2 + 4k^2}}{j + l + \sqrt{(j - l)^2 + 4k^2}}}
$$

$$
\begin{align}
j &= \frac{1}{2}\left[\frac{1}{q^2_1 + q^2_2}\right] - \frac{\sin{\theta}^2\cos{\psi}^2(q^2_1 + q^2_2 - 2)}{q^2_1 q^2_2} \\& - \left( \frac{1}{q_1^2}- \frac{1}{q_2^2}\right)(\cos{2\varphi}(\cos^2{\theta}\cos{\phi}^2 - \sin^2{\psi}) - \cos{\theta}\sin{2\varphi}\sin{2\psi})
\end{align}
$$
$$
\begin{align}
k &= \frac{1}{4 q_1^2 q_2^2}\Big[2\cos{\theta}(q_1^2 - q_2^2)\cos{2\phi}\cos{2\varphi} + \\&\Big(\sin{\theta}^2(q_1^2 + q_2^2 - 2) + (1 + \cos{\theta}^2)(q_1^2 - q_2^2)\cos{2\varphi} \Big)\sin{2\psi}\Big]
\end{align}
$$
$$
\begin{align}
l &= \frac{1}{2}\Big[\left(\frac{1}{q_1^2} + \frac{1}{q_2^2}\right) - \frac{\sin{\theta}^2\sin{\phi}^2(q_1^2 + q_2^2 - 2)}{q_1^2 q_2^2} - \\ & \left(\frac{1}{q_1^2} - \frac{1}{q_2^2}\right)\Big(\cos{2\varphi}(\cos{\theta}^2 \sin{\phi}^2 - \cos{\psi}^2) + \cos{\theta}\sin{2\varphi}\sin{2\phi}\Big)\Big]
\end{align}
$$

The orientation angle in the plane of the sky is
$$
\theta_{\epsilon} = \tan^{-1}{\left[\frac{l - j + \sqrt{(j - l )^2 + 4k^2}}{2k}\right]}
$$
And the elongation parameter
$$
\epsilon_{\parallel} = \sqrt{\frac{q_p}{q_1q_2}}f^{-3/4}
$$
$$
f = \left[\sin{\theta}^2\left(\frac{\sin{\varphi}^2}{q_1^2}+\frac{\cos{\varphi}^2}{q_2^2}\right) + \cos{\theta}^2\right]
$$
the semimajor axis of the projected ellipse is 
$$
l_p = \frac{l_s}{\epsilon_{\parallel}\sqrt{f}}
$$
and the projected lenght scales $l_s$ and $l_{\mathrm{los}}$ are related by
$$
l_{\mathrm{los}} = l_s/\sqrt{f}
$$
Finally the elliptical radius on the plane of the sky becomes
$$
\xi^2 = \left(x^{2}_1 + \frac{x^2_2}{q_p^2}\right)\left(\frac{l_s}{l_p}\right)^2
$$
The measured 2-D projected profiles is given by
$$
F_{2D}(x_{\xi},l_p |p) = 2l_p\epsilon_{\parallel}\int_{x_{\xi}}F_{\mathrm{3D}}(x_{\zeta}, l_s| p)\frac{x_{\xi}}{\sqrt{x_{\zeta}^2 - x_{\xi}^2}}dx_{\zeta}
$$
Where $x_{\zeta} = \zeta/l_s$ and $x_{\xi} = \xi/l_p$. $l_s$ is obtained by eigvalues of shape matrix
$$
T = \text{diag}{[1, q_1, q_2]}
$$

