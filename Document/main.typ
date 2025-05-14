// #import "@preview/light-report-uia:0.1.0": *
#import "lib.typ": report
// CHANGE THESE
#show: report.with(
  title: "Properties of Pair of Coupled Maps",
  authors: (
    "Erin Jossy",
  ),
  course_code: "AM5650",
  course_name: "Non Linear Dynamics",
  date: "August 2024",
  lang: "en", // use "no" for norwegian
)

// neat code
#import "@preview/equate:0.3.1": equate

#show: equate.with(breakable: true, number-mode:"label")
#set math.equation(numbering: "(1)",supplement: [Eq.])

#import "@preview/equate:0.3.1"

= Introduction
The natural world and engineered systems are filled with phenomena governed by nonlinear interactions. Unlike linear systems where outputs are proportional to inputs and superposition holds, nonlinear systems exhibit a rich tapestry of behaviors, including multiple stable states, periodic oscillations, quasi-periodicity, and the hallmark of deterministic chaos – extreme sensitivity to initial conditions alongside bounded, aperiodic long-term evolution. Understanding and predicting the behavior of these nonlinear systems is a central challenge across diverse scientific and engineering disciplines. Coupled map systems consist of networks of relatively simple dynamical units (maps) that interact according to specific coupling rules, giving rise to remarkably rich collective behaviors. Originally introduced as mathematical abstractions, coupled maps have evolved into powerful modeling tools with applications spanning physics, biology, neuroscience, economics, and numerous other fields.

In reality, few systems exist in complete isolation. More often, individual units or subsystems interact with each other, forming coupled systems. The study of coupled dynamical systems investigates how the behavior of individual components is modified by their mutual influence, and critically, what new collective or emergent behaviors arise from these interactions. This coupling can take many forms, from diffusive coupling representing the exchange of mass or energy, to all-to-all coupling in certain network structures, or more complex, specific interaction topologies.

What makes coupled map systems particularly intriguing is their ability to generate complex global behaviors from simple local dynamics and interactions. A single map, such as the logistic map $x_(n+1) = r x_n (1-x_n)$, can already exhibit period-doubling bifurcations leading to chaos. When multiple such maps are coupled together, new phenomena emerge: synchronization, where initially different units begin to behave identically; pattern formation, where spatial structures spontaneously develop; and various forms of intermittency and turbulence that can model real-world phenomena across scales.
$ x_(n+1) = f(x) $
If an initial condition $x_0$ is specified at time level n = 0, the system state at time level n = 1 is $x_l = f(x_0)$ ,and so on. The function f can be any linear or non-linear map, for instance, the exponential map, $ f(x) =x e^(r(1 - x ) ) #<function> $  as per Ricker to model the population dynamics of a prey species
In this paper, explores the properties of a pair of coupled logistic maps given by the following equations:

$ x_(n+1) &= d f(x_n) + (1-d) f(y_n) \  
 y_(n+1) &= (1-d) f(x_n) + d f(y_n) #<equation> $

where $d$ is the coupling strength, which is a number between 0 and 1 and represents the fraction of the total population that interacts with the other population. This form of coupling often arises in physical systems, such as in the study of coupled oscillators. The equation can be interpreted as follows: at each time step $n$, the population $x_n$ is updated by taking a weighted average of its own value and the value of the other population $y_n$ with weights $d$ and $(1-d)$ respectively. The weights are chosen such that the total population remains constant. For example, if $d=0.5$, the population $x_n$ is updated by taking the average of its own value and the value of $y_n$. This form of coupling is often used to model the interactions between two species, such as predator-prey systems or competing species. The authors investigate several "global, generic properties" that emerge from this coupling structure, independent of the specific functional form of f (though their numerical examples use specific maps).
We shall denote the mapping  for
brevity, as
$  (x_(n+1),y_(n+1)) = M(d) circle.stroked.small (x_n,y_n) $

The paper have revealed that coupled maps has a great number of phenomena, including synchronization, chaos, and the emergence of complex patterns. These systems can exhibit a wide range of behaviors depending on the coupling strength, the nature of the maps, and the initial conditions. For example, weak coupling can lead to synchronization, where the maps converge to a common trajectory, while strong coupling can induce chaotic behavior or even lead to the formation of stable patterns.

This report aims to provide a comprehensive examination of coupled map systems, beginning with fundamental theoretical concepts and progressing through dynamical properties to diverse applications. We will explore  established results, highlighting the universal principles that govern these systems while acknowledging the domain-specific insights they have yielded. Throughout, we emphasize the mathematical elegance of coupled maps and their practical utility as models of complex, real-world phenomena. 


The investigation of these global properties is significant because it provides a framework for understanding how interaction shapes collective behavior. Discovering that simple coupling can lead to:
- predictable symmetries (reducing analytical effort),
- the emergence of stable periodic states from chaotic components (offering mechanisms for control or explaining observed stability in nature), and
- robust synchronized states (with implications for communication and coordinated activity)
The paper's @UDWADIA199816 findings highlight that even in seemingly simple two-unit coupled systems, a rich and structured array of dynamical behaviors can emerge, driven by the interplay between the intrinsic dynamics of the individual units and the nature of their coupling. The "seven-zone dynamics" characterized by distinct bifurcations (tangent, Sacker/Hopf-like) further underscores the structured complexity arising from these interactions.

= Analytical Results
Some analytical results are discussed below, which will be useful in understanding the properties of coupled maps. The results are derived from the equations @equation  the coupled maps and are presented in a concise manner.
1. *Result 1*: For each orbit ${(hat(x)_n, hat(y)_n) | n = 0,1,2 . . . . }$ corresponding to the parameter value d, there corresponds an orbit ${(hat(y)_n, hat(x)_n) | n = 0, 1,2 . . . . }.$

  *Corollary 1*: For each n-periodic orbit of the map
  described by ${(hat(x)_n, hat(y)_n) | n = 0,1,2 . . . .}$ corresponding to the parameter d, there exists another $n$-periodic
  orbit described by ${(hat(y)_n, hat(x)_n) | n = 0, 1,2 . . . . }.$

2. *Result 2*: Consider the orbit ${(hat(x)_n,hat(y)_n) | n = 0, 1, 2,
...}$ corresponding to a certain value of the parameter $d = ½- d_o, 0 < d_o < ½$. For each such orbit, there corresponds an orbit given by ${(hat(x_0), hat(y_0)),(hat(y)_1, hat(x)_1),(hat(x)_2, hat(y)_2),(hat(y)_3, hat(x)_3) . . . . }$, corresponding to the parameter $d = ½ + d_0$, i.e., each alternate point of the second orbit has its x- and y-coordinate switched with respect to the corresponding point of the first orbit.
  
  *Corollary 2*: If the map $M(½ - d_0), 0 <= d_0 <= ½$, has a $2n$-period orbit starting from some $(x_0, y_0), n =
  1,2 . . . . .$ then the map $M(½ + d_0)$ must have a $2n$-
  period orbit starting from the same point $(x_0, y_0), n =
  1,2 . . . .$

  *Corollary 3*: If the map $M(½ - d_0), 0 <= d_0 <= ½$, has
  a $2n-1$-period orbit starting from some $(x_0, y_0), n =
  1,2 . . . . .$ with $x_0 != y_0$ then the map $M(½ + d_0)$ must have a $2(2n-1)$-
  period orbit starting from the same point $(x_0, y_0), n =
  1,2 . . . .$

3. *Result 3*:The orbit of the map starting with $(x_0, x_0)$ will always consist of points of the form $(x_n, x_n)$. In other words, orbits which begin on the diagonal in the $(x, y)$ phase space lie entirely on the diagonal.
4. *Result 4*: For $d = ½$ the orbit of the map starting with $(x_0, y_0)$ will consist of points of the form $(x_n, x_n)$ after the first iteration.
5. *Result 5*: Each Lyapunov exponent of the map $M(½ -d_0), 0 <= d_0 <= ½$, starting from some $(x_0, y_0)$ is the same as that of the map $M(½ +d_0)$ starting from the same point $(x_0, y_0)$.
6. *Result 6*: When $d = l$, one Lyapunov exponent tends and to $-infinity$.

These analytical results provide a robust framework for understanding the behavior of the coupled map system. They reveal inherent symmetries related to variable interchange and the coupling parameter d, explain the special invariant nature of the diagonal subspace, connect the coupled system's dynamics directly to the underlying 1D map f along the diagonal and at $d=0.5$, demonstrate the mirrored stability properties around $d=0.5$ via Lyapunov exponents, and analytically predict the boundaries of the stable synchronization regime. These findings hold generally for the given coupling structure and significantly simplify the interpretation of complex numerical observations.

= Numerical Results
The numerical results are obtained by simulating the coupled map system @equation using the equations governing the dynamics. The simulations are performed for different values of the coupling strength d and the initial conditions. The results are presented in the form of plots and tables, which illustrate the behavior of the system under different conditions.

@map show the the behavior of the system for different values of the coupling strength d. The bifurcation diagram shows the periodic behavior of the system for certain ranges of d, while the Lyapunov exponent plot @lyp shows the stability of the system for different initial conditions. The results indicate that the coupled map system exhibits a rich variety of behaviors, including periodic orbits, chaotic behavior, and synchronization. Thus coupling two
chaotic systems can stabilize both of them. In fact,
as @map shows, for 0.03 < d < 0.13 we observe two
one-period orbits, i.e., by coupling two chaotic units
we can arrive at complete equilibrium. This
stabilizing phenomenon may very well explain why
there is so much stability in this physical world despite
all the reported chaos.

#show figure.caption: (t)=>align(left,t)

#figure(
  grid(
    rows: 2,
    image("/Graphics/1_a.jpg", width: 40%),
    image("/Graphics/1_b.jpg", width: 80%),

  ),
    caption: [ _Bifurcation diagram of the exponential map plotted against the growth parameter $r$. For r say, r = 4, the map is chaotic. By coupling two such maps in the form described by equation, we can obtain periodic behavior for certain ranges of the coupling parameter d._],
)<map>
When synchronization is observed, the response of the coupled system, after a large number of iterations, is limited to
the line x = y. Thus, instead of spanning the two- dimensional surface, the dynamics is confined to a one-dimensional line. This one-dimensional dynamics may be viewed as a step towards stabilization.

When the two units evolve synchronously, both populations are identical, i.e., $x = y$. Therefore, the difference between the two populations $x - y$ should be identically zero. Synchronicity variation with the coupling parameter d, could be understood from @sync.@Raju1997 It can be seen that the range of d for which the dynamics is synchronous is significant.
coupled exponential map displays a variety of behavior over the range of d values of interest. The overall dynamics of the coupled exponential map
system can broadly be divided into seven zones as in @sync(c). The transition from one zone to another is marked by a bifurcation with respect to the parameter d. Either a tangent bifurcation(Zone II to Zone III) or a Hopf-like bifurcation (Zone II to Zone III) may mark the transition from one zone to another. 

Some features of dynamics in the seven zones are as follows:
- Zone I: The system is primarily chaotic and desynchronized. The two maps follow different chaotic paths, and the largest Lyapunov exponent is positive.
- Zone II: The system is stable and periodic, with the two maps converging to a finite set of points. Both Lyapunov exponents are negative, confirming stability.
- Zone III: The system is chaotic and desynchronized, with trajectories filling a complex region similar to Zone I. The largest Lyapunov exponent is positive.
- Zone IV: The system is chaotic and synchronized, with trajectories starting off the diagonal x=y being attracted to it. The tangential Lyapunov exponent is positive, while the transverse Lyapunov exponent is negative.
- Zone V: The system is chaotic and desynchronized, similar to Zone III. The largest Lyapunov exponent is positive.
- Zone VI: The system is stable and periodic, with trajectories converging to a finite set of points. Both Lyapunov exponents are negative.
- Zone VII: The system is primarily chaotic and desynchronized, similar to Zone I. The largest Lyapunov exponent is generally positive.

#figure(
  grid(
    rows: 2,
    image("/Graphics/2_a.jpg", width: 90%),
    image("/Graphics/2_d.jpg", width: 50%),
  ),
  
    caption: [ _Plots depicting the diagonal attraction  A coupled exponential map with $r = 4$ is used: $(a) d = 0.16$; $(b) d = 0.21$. For some values of $d$, as in (a), the
    response is over a region in the line $x = y$. (c) The difference $x - y$ plotted against $d$. Over a large range of $d$, the coupled response
    is synchronous as seen by the null values. (d) The transitions from zone 2 to zone 1, and from zone 6 to zone 7 are marked by Hopf bifurcations $(r = 4)$. The closed loop trajectory on the left is caused by a bifurcation of the one-period orbit in zone 2. The corresponding closed trajectory, caused by a Hopf bifuraction of the two-period orbit in zone 6, is shown on
    the right. This orbit alternates between the two closed loops. The loops are symmetric about the line x = y.._],
)<sync>
One of the major question would be what is the probability of the two maps being synchronized. The probability of synchronization is defined as the fraction of time that the two maps are synchronized. The probability of periodicity, synchronicity, and their total are plotted in @prob for cases where the underlying single map is chaotic.
#figure(
    image("/Graphics/prob.jpg", width: 40%),
  
    caption: [_The probability of obtaining (i) periodic behavior (dashed lines), (ii) synchronous behavior (dashed dot) and (iii) their total,
    plotted against $r$._],
)<prob>
#figure(
    image("/Graphics/Lyp.jpg", width: 60%),
  
    caption: [_The spectrum of Lyapunov exponents. (a) The ICs do not fall on the diagonal $(x_0!=y_0)$
    here are ranges of $d$ for which both exponents are negative which indicates that the behavior in those ranges is non-chaotic. (b) The
    ICs fall on the diagonal $(x_0 = y_0)$. At least one exponent is  positive indicating that the diagonal is chaotic._],
)<lyp>
= Novel results
== Largest Lyapunov Exponent map
2D parameter scan of Lyapunov exponents for the coupled exponential map. The largest Lyapunov exponent is plotted against the coupling strength d and map non-linearity r .Lines or sharp transitions in the color map corresponding to $λ₁ ≈ 0$ indicate bifurcation curves in the $(r, d)$ parameter plane.@ly_map helps understand how the seven zones discussed in the paper for a fixed r (like r=4) evolve or change shape as r itself is varied.
#figure(
    image("/Graphics/lyp_map.jpg", width: 50%),
  
    caption: [_Contour map of the largest Lyapunov exponent for a coupled exponential map._],
)<ly_map>
== Synchronization Time vs. Coupling d
@time, Shows how quickly the system synchronizes once the coupling d is within the stable synchronization zone (Zone IV). offers valuable insights into the dynamics and efficiency of the synchronization process itself. If rapid synchronization is a desired feature (e.g., in secure communication systems where quick lock-in is needed, or in neural systems for fast processing), this plot can help identify the range of coupling strengths d that minimize the time to synchronize. 


#figure(
    image("/Graphics/time.jpg", width: 50%),
  
    caption: [_Show the iterations required for synchronization for various values of d._],
)<time>
= Conclusion
Analytical derivations and numerical simulations, we confirmed the paper's findings on quasi-symmetry around d=0.5, the emergence of periodicity, and the significant phenomenon of chaotic synchronization within distinct zones of coupling strength. Lyapunov exponent calculations validated the stability and chaotic nature of these zones. Furthermore, additional visualizations, such as 2D parameter scans of the largest Lyapunov exponent and synchronization time analyses, were proposed and demonstrated to offer deeper global insights into the system's rich dynamics, highlighting the profound impact of coupling strength and intrinsic map nonlinearity on the collective behavior of even simple interacting systems. These results underscore the power of coupling to induce order or complex synchronized states from potentially chaotic components, with broad implications for understanding natural and engineered systems.
// == Tables
// Here's a table:
// #figure(
//   caption: [Table of numbers],
//   table(
//     columns: 6,
//     inset: 10pt,
//     align: horizon,
//     table.header(
//       [*Zone*], [*Dynamics*],
//       [*Phase Space*], [*Lyapunov Exponents*],
//       [*x-y Plot*], [*Bifurcation*] 
//     ),
//     [Zone I ($0< d< 0.03$ for r=4)], [ Primarily Chaotic, Desynchronized.],[$x_n$and $y_n$ follow different chaotic paths.],
//     [Largest LE is +ve, transverse LE is also typically +ve. ], [Shows a broad band of non-zero values.],[Sacker (Hopf-like) bifurcation occurs at boundary from Zone II to Zone I.],
//     [Zone II ($0.03 ≤ d ≤ 0.13$ for r=4)], [Stable Periodic, Desynchronized.], [Trajectories converge to a finite set of points.],
//     [Both Lyapunov exponents are -ve, confirming stability.], [Settles to one of two distinct, non-zero constant values.],[Tangent bifurcation occurs, stable fixed points are annihilated.],
//     [Zone III ($0.13 < d < 0.21$ for r=4)], [Chaotic, Desynchronized.], [Trajectories again fill a complex region, similar to Zone I but potentially with a different structure.],
//     [Largest LE is +ve.], [Shows a broad band of non-zero, fluctuating values.],[Transverse LE crosses zero from +ve to -ve.],
//     [Zone IV ($0.21 ≤ d ≤ 0.79$ for r=4)], [Chaotic, Synchronized.], [Trajectories starting off the diagonal $x=y$ are attracted to it.],
//     [Tangential LE is +ve, transverse LE is -ve.], [Quickly collapses to zero and stays there.],[Transverse LE crosses zero from -ve to +ve.],
//     [Zone V ($0.79 < d < 0.87$ for r=4)], [Chaotic, Desynchronized.], [Similar to Zone III.],
//     [Largest LE is +ve.], [Shows a broad band of non-zero, fluctuating values.],[Re-emergence of stable periodic behavior, likely through inverse bifurcations.],
//     [Zone VI ($0.87 ≤ d ≤ 0.97$ for r=4)], [Stable Periodic, Desynchronized.], [Trajectories converge to a finite set of points.],
//     [Both Lyapunov exponents are -ve.], [Settles to a pattern that alternates between two distinct non-zero values.],[Sacker (Hopf-like) bifurcation occurs, stable period-2 orbit loses stability.],
//     [Zone VII ($0.97 < d < 1$ for r=4)], [Primarily Chaotic, Desynchronized.], [Similar to Zone I.],
//     [Largest LE is generally +ve.], [Shows a broad band of non-zero values.],[This zone extends to d=1 (uncoupled, independent chaotic maps).],
//   ) 
// )
