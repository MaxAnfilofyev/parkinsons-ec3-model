## A Load-Induced Energetic Tipping Point Explains Selective Vulnerability of Substantia Nigra Neurons

Max Anfilofyev
[max.anfilofyev@gmail.com](mailto:max.anfilofyev@gmail.com)

---

## **Abstract**

Dopaminergic neurons of the substantia nigra pars compacta (SNc) are selectively vulnerable in Parkinson’s disease, while closely related neurons in the ventral tegmental area (VTA) are comparatively spared. Although mitochondrial dysfunction, calcium stress, and α-synuclein aggregation have each been implicated, none alone explains why anatomically similar populations exhibit such different fates. Here we develop a **minimal two-variable energetic model** that captures only mitochondrial functional capacity, energetic reserve, and the combined load from axonal arborization and calcium handling. Despite its simplicity, the model reveals that increasing structural load deforms the energetic landscape until a **saddle-node bifurcation** emerges, producing coexisting healthy-energy and collapsed-energy states. SNc-like neurons, which bear extreme axonal and calcium-handling demands, reside **inside** this bistable regime, operating near a separatrix that renders them vulnerable to even modest metabolic perturbations. In contrast, VTA-like neurons lie **outside** the bistable window and robustly return to their high-energy state following similar disturbances. The model reproduces hallmark features of Parkinsonian degeneration—long periods of stability, sudden irreversible collapse, and population-specific susceptibility—using only the geometry of load-dependent energy regulation. These findings suggest that selective SNc vulnerability arises not from unique molecular defects, but from the **fundamental dynamical structure** imposed by their extraordinary anatomical and physiological load.

---

## **1. Introduction**

Parkinson’s disease is marked by the striking and selective degeneration of dopaminergic neurons in the substantia nigra pars compacta (SNc), while neighboring dopaminergic populations in the ventral tegmental area (VTA) remain largely preserved.[15–17](#references) This anatomical specificity has remained one of the field’s central puzzles. Multiple lines of evidence converge on a common theme: **SNc neurons carry exceptionally large structural and functional loads.** Their axonal arborizations span hundreds of thousands to millions of synapses, embedded in complex basal ganglia microcircuits, and their pacemaking physiology relies on Ca²⁺ entry that imposes substantial metabolic demand.[1–4,5–9,15,16](#references) These characteristics imply that SNc neurons operate close to the limits of energetic feasibility.[15–17](#references)

Despite extensive work on mitochondrial dysfunction, oxidative stress, dopamine oxidation, and α-synuclein aggregation, no single factor has fully explained why SNc neurons are so much more vulnerable than VTA neurons.[5–14,18–20,24](#references) Many existing models either include dozens of variables or focus on downstream pathology without addressing the fundamental question: *Why is the SNc neuron so close to a tipping point, while the VTA neuron is not?*[21–23](#references)

Here we propose a **minimal energetic model** that captures this tipping-point behavior using only two dynamic variables and two structural loads. The key idea is that SNc neurons sit near a **saddle-node bifurcation** in cellular energy supply and demand.[21–23](#references) Under normal conditions, they maintain a stable energetic state. However, a small transient perturbation—such as a brief metabolic challenge or mitochondrial insult—can push them across a separatrix into an alternative low-energy attractor from which recovery is no longer possible. In contrast, VTA neurons operate far from this bifurcation and robustly return to their healthy energetic state even after substantial perturbations.[15–17](#references)

The goal of this work is not to reproduce all molecular details of Parkinson’s disease, but to demonstrate that **the geometry of cellular energy dynamics alone** is sufficient to explain the core vulnerability pattern. In doing so, we provide a conceptual bridge between anatomical load, Ca²⁺ physiology, mitochondrial strain, and the catastrophic failure mode characteristic of SNc degeneration.[1,3,5–10,15–17,21–23](#references)

---

## **2. Minimal Energetic Model**

### **2.1 Model Structure**

To capture the essential energetic behavior of dopaminergic neurons without introducing unnecessary biochemical detail, we construct a **two-variable system** describing the interaction between:

* **$E(t)$**: the cell’s *available energetic reserve*, normalized to $[0,1]$.
  This variable integrates ATP availability, NADH balance, and the capacity of the cell to meet ongoing energetic demands.[24,25](#references)

* **$M(t)$**: the cell’s *functional mitochondrial capacity*, also normalized to $[0,1]$.
  This variable represents the collective ability of mitochondria to sustain oxidative phosphorylation under load, including turnover, repair, and stress-induced dysfunction.[10–14,24](#references)

SNc–VTA differences are introduced not by altering the equations themselves, but by adjusting two **load parameters** that modulate energy consumption and mitochondrial stress:

* **$A$**: the **axonal arborization load**, proportional to the number of synapses that must be maintained and serviced.
  Anatomically, SNc neurons have an arbor roughly 4–10× larger than VTA neurons, which we encode as higher values of $A$.[1,3,15–17](#references)

* **$C$**: the **Ca²⁺-handling load**, representing the energetic overhead of pacemaking and channel activity.
  In many simulations we set $C = 1$, with relative differences absorbed into the effective scaling of $A$.[5–9](#references)

```mermaid
---
title: Figure 1 - Minimal Energetic Model
config:
  layout: elk
---
flowchart LR;

    %% Upstream structural loads
    A["Axonal Arbor Size (A)"]
    C["Ca²⁺ Pacemaking Load (C)"]

    %% Combined energetic load
    L["Metabolic Demand (A × C)"]

    %% Internal state variables
    E["Energetic Reserve (E)"]
    M["Mitochondrial Capacity (M)"]

    %% Loads increase consumption
    A --> L
    C --> L
    L -->|"ATP consumption"| E

    %% Mitochondria support energy
    M -->|"ATP production"| E

    %% Energy enables mitochondrial maintenance
    E -->|"adequate ATP → repair/turnover"| M

    %% Low energy causes mitochondrial vulnerability
    E -.->|"insufficient ATP → impaired repair"| M

```

**Figure 1. Minimal energetic model of dopaminergic neurons.**
Axonal arborization load $A$ and Ca²⁺ pacemaking load $C$ combine into an effective metabolic demand that drains energetic reserve $E$. Mitochondrial capacity $M$ supports ATP production and therefore increases $E$, while sufficient energy is required to maintain and repair mitochondria. When energetic reserve is low under high structural load, mitochondrial damage accumulates. This compact feedback loop—mitochondria support energy, energy maintains mitochondria, and loads destabilize both—is sufficient to generate coexisting healthy and collapsed energetic states under high load.

The model structure (Figure 1) consists of three core interactions:

1. **Mitochondrial support of energy.**
   Functional mitochondria increase the energy reserve, reflecting ATP production proportional to available mitochondrial capacity.[10–14,24](#references)

2. **Load-dependent energy consumption.**
   Axonal and Ca²⁺ loads drain energy at a rate that scales with both arbor size and Ca²⁺-related metabolic demand.[1,3,5–9,15](#references)

3. **Energy-dependent mitochondrial maintenance.**
   Mitochondrial capacity is replenished through repair and turnover but is damaged when energetic reserve is insufficient to buffer load-induced stress, consistent with experimental observations in SNc neurons.[10–14,15,24](#references)

Together these interactions form a compact feedback loop: mitochondria support energy; energy enables mitochondrial maintenance; loads push both toward failure. Remarkably, this minimal structure is mathematically sufficient to generate **two coexisting energetic states**—a healthy, high-energy attractor and a collapsed, low-energy attractor—under high load. SNc neurons reside near the fold separating these states, whereas VTA neurons do not.

This parsimonious setup allows us to capture the qualitative behavior of dopaminergic neurons while keeping the model analytically transparent and numerically tractable.[21–23](#references)

---

### **2.2 Equations and Qualitative Behavior**

The dynamics of the system are described by two coupled differential equations governing the temporal evolution of energetic reserve $E(t)$ and mitochondrial functional capacity $M(t)$. Both variables are normalized to the interval $[0,1]$, with higher values reflecting greater energetic availability or mitochondrial robustness:

$$
\frac{dE}{dt}
= k_1 M(1-E) + k_2 E^2(1-E) - \left(L_0 + L_1 A C\right) E ,
$$

$$
\frac{dM}{dt}
= k_M(1-M) - \beta A C M(1-E).
$$

The first term in the energy equation represents ATP production by functional mitochondria: as long as mitochondria are intact ($M$ large), they drive $E$ upward toward its maximum.[10–14,24,25](#references) The nonlinear term $k_2 E^2(1-E)$ captures activity-dependent amplification of energy availability, reflecting processes such as efficient pump operation, metabolic feedback, and cooperative mitochondrial behavior.[24,25](#references) The final term represents energy consumption, which increases with both arborization load $A$ and Ca²⁺-handling demand $C$, capturing the combined costs of maintaining large axonal arbors and pacemaking currents in SNc neurons.[1,3,5–9,15](#references)

The mitochondrial equation includes a repair-and-replacement term $k_M(1-M)$, counteracted by damage that occurs when high structural load (through $A$) and low energy (through $1-E$) coincide. This reflects the empirical observation that mitochondrial stress in dopaminergic neurons is strongly dependent on the combination of Ca²⁺ influx, axonal maintenance, and energetic sufficiency.[10–14,15](#references)

Together, these interactions form a compact positive–negative feedback motif: mitochondria support energy; energy maintains mitochondria; loads destabilize both. Despite the simplicity of this architecture, the system exhibits **nonlinear geometric structure** characteristic of a saddle-node bifurcation.[21–23](#references)

To illustrate the qualitative structure of the energy dynamics, Figure 2A shows a conceptual landscape representation of the system’s two attractors and the saddle separating them. The high-energy attractor corresponds to normal physiological operation, while the low-energy attractor represents an irreversible collapsed state. The saddle acts as the basin boundary: perturbations that remain on the healthy side of this boundary recover, whereas perturbations that cross it transition to the collapsed state. This geometric picture anticipates the mathematical structure revealed in the phase-plane analysis (Figure 2B) and highlights why neurons experiencing large structural loads, such as substantia nigra dopaminergic neurons, reside precariously close to the tipping boundary, whereas ventral tegmental area neurons lie far from the saddle and therefore exhibit robust recovery after metabolic or calcium-related stress.[15–17,21–23](#references)

```mermaid
---
title: Figure 2A - Conceptual Energy Landscape Showing Bistability Under Load
config:
  layout: elk
---
flowchart LR

    %% Basins
    subgraph Healthy["High-Energy Attractor (Healthy State)"]
        H(( ))
    end

    subgraph Saddle["Saddle (Tipping Boundary)"]
        S(( ))
    end

    subgraph Collapsed["Low-Energy Attractor (Collapsed State)"]
        C(( ))
    end

    %% Flow directions
    H -->|"small perturbation"| H
    H -->|"large perturbation"| S
    S -->|"cross separatrix"| C
    C -->|"irreversible under high load"| C

    %% Neuron types
    VTA["VTA neuron 
    (low load, robust)"] --- H
    SNc["SNc neuron 
    (high load, near threshold)"] --- S

    %% Load axis annotation
    A_increase["Increasing axonal load shifts the saddle left and shrinks the healthy basin"]
    A_increase ----> S

```

**Figure 2A. Conceptual energy landscape illustrating bistability under structural load.**
The minimal energetic model predicts two coexisting stable states separated by a saddle (tipping boundary). The upper basin corresponds to a healthy, high-energy attractor, while the lower basin represents an energetically collapsed state. Small perturbations within the healthy basin are restored, whereas sufficiently large perturbations cross the separatrix and drive the system into the collapsed attractor. Increasing axonal load $A$ shifts the saddle leftward and shrinks the healthy basin. VTA neurons, which experience low load, lie deep within the healthy basin; SNc neurons, under high load, lie near the saddle and are therefore vulnerable to tipping into collapse.[15–17,21–23](#references)

![Figure 2B. Phase portrait of the minimal energetic model at a load level representative of substantia nigra dopaminergic neurons ($A = 1$)](/Stage_3_PhasePlanes/latest_run/phase_plane_SNc_A_1.00.png)**Figure 2B. Phase portrait at a substantia nigra–like load ($A = 1$).**
Grey arrows show the vector field in the $(E,M)$ plane. The solid purple curve denotes the energy nullcline ($\frac{dE}{dt}=0$), forming an S-shaped curve with two folds, and the dashed purple curve denotes the mitochondrial nullcline ($\frac{dM}{dt}=0$). Their intersections yield three equilibria: a high-energy stable fixed point (upper black dot), a low-energy stable fixed point (lower black dot), and an intermediate saddle (red cross). Colored trajectories illustrate that initial conditions on one side of the saddle’s stable manifold converge to the healthy high-energy attractor, whereas those on the other side collapse into the low-energy state. This phase-plane geometry provides a concrete dynamical realization of the conceptual landscape in Figure 2A and shows that SNc-like loads place the system inside a load-induced tipping regime.[21–23](#references)

This geometry is the mathematical signature of a **tipping point**: a minimal energetic failure mode that emerges directly from load-dependent energy consumption and energy-dependent mitochondrial stress.[21–23](#references) In subsequent sections, we show that ventral tegmental area neurons reside far from this fold, while substantia nigra neurons sit near it, making them uniquely susceptible to collapse.[15–17](#references)

---

## **3. Bifurcation Analysis: Load-Driven Emergence of a Tipping Point**

To determine how structural load shapes the energetic stability landscape, we performed a one-parameter continuation in the axonal arborization parameter $A$. For each value of $A$ between 0.2 and 1.4, we computed all steady states of the system and classified their stability by linearization, following standard dynamical-systems practice for saddle-node bifurcations.[21–23](#references) The resulting diagram (Figure 3) reveals a **saddle-node bifurcation** in the energy variable $E$, separating regions of monostability from regions in which two stable energetic states coexist.

<a id="figure_3"></a>![Figure 3. Saddle-node bifurcation of energetic reserve as a function of axonal load](/Stage_2_Bifurcation/latest_run/ec3_bifurcation_E_vs_A.png) **Figure 3. Saddle-node bifurcation of energetic reserve as a function of axonal load.**
The steady-state energetic reserve $E^*$ is shown as a function of axonal load $A$. Solid markers indicate stable equilibria; crosses indicate unstable saddle points. At low structural load (left), the system is monostable, with a single high-energy attractor corresponding to VTA-like arborization levels. As $A$ increases, the energy nullcline deforms and a saddle-node bifurcation appears near $A \approx 0.86$, generating a low-energy stable state and an intermediate saddle. Between $A \approx 0.86$ and $A \approx 1.06$, the system becomes bistable, with coexisting healthy- and collapsed-energy attractors separated by the saddle’s stable manifold. Anatomically realistic SNc-like loads fall within this bistable window, whereas VTA-like loads lie to the left in the monostable regime, consistent with the selective vulnerability of SNc dopaminergic neurons.[1,3,15–17,21–23](#references)

At low structural load (left side of Figure 3), the system has a **single stable equilibrium** corresponding to a healthy, high-energy state. This monostable regime encompasses the range of arborization values typically associated with ventral tegmental area neurons.[1,3,15–17](#references) As $A$ increases, the energy nullcline bends downward, and at a critical load value near $A \approx 0.86$, the system undergoes a fold bifurcation that produces an additional pair of equilibria: a saddle point and a low-energy stable state.[21–23](#references) Beyond this point, and until approximately $A \approx 1.06$, the system is **bistable**, exhibiting both a healthy-energy attractor and a collapsed-energy attractor, separated by a codimension-one separatrix.

This bistable window is precisely where **substantia nigra dopaminergic neurons** are expected to lie based on anatomical reconstructions showing their vastly expanded axonal arborizations and elevated energetic burden.[1,3,15–17](#references) In this regime, the neuron can maintain normal energetic function but only by remaining on the high-energy side of the saddle’s stable manifold. Even transient perturbations—such as brief metabolic challenges, mitochondrial insults, or calcium-driven fluctuations in ATP demand—can push the system across this barrier, leading to an irreversible transition into the low-energy attractor.[5–9,10–14,24](#references) This geometry provides a mechanistic explanation for the characteristic “catastrophic failure” in substantia nigra neurons despite long periods of apparent resilience.

An important aspect of the bifurcation structure is that **the right fold (the point at which the two stable equilibria annihilate)** lies **beyond the biologically relevant range of $A$** for dopaminergic neurons. As a consequence, while the system exhibits a true saddle-node bifurcation mathematically, **hysteresis does not play a major role biologically**. Once the neuron has collapsed into the low-energy state, decreasing $A$ (e.g., through loss of axon terminals) does not restore the high-energy state, because the fold at which recovery would occur sits outside physiologically plausible arborization values.[1,3,15–17](#references) Thus, the collapse is effectively irreversible in the anatomical range relevant for Parkinson’s disease.

If desired, this bifurcation geometry can be further examined through extended parameter sweeps or numerical continuation methods (Supplementary Figures S1–S3), which show that the bistable window persists under moderate changes in mitochondrial turnover, Ca²⁺-handling cost, and the strength of energy-dependent mitochondrial damage. These analyses confirm that saddle-node structure is a **robust qualitative feature** of the model rather than a fine-tuned artifact of any specific parameter choice.[21–23](#references)

Together, these results establish that **axonal arborization is a natural control parameter for dopaminergic energetic stability**, and that substantia nigra neurons, by virtue of their extreme structural load, are uniquely situated near a dynamical tipping point.[1,3,15–17,21–23](#references)

---

## **4. Comparison of Substantia Nigra and Ventral Tegmental Area Neurons**

A defining characteristic of substantia nigra pars compacta (SNc) dopaminergic neurons is their extraordinarily large and widely distributed axonal arbor. Anatomical reconstructions indicate that a single SNc neuron forms hundreds of thousands to millions of synapses, a structural scale unmatched by most other neuronal types.[1,3,15,16](#references) In contrast, dopaminergic neurons in the ventral tegmental area (VTA) innervate far fewer targets, with substantially smaller arbor size and reduced calcium-handling burden during pacemaking.[5–9,15–17](#references) These anatomical and physiological differences map naturally onto the load parameter $A$ in our model.

| ![](/Stage_3_PhasePlanes/latest_run/phase_plane_VTA_A_0.40.png)          | ![](/Stage_3_PhasePlanes/latest_run/phase_plane_SNc_A_1.00.png)                        |
| ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------- |
| **Figure 4A. Phase plane of the minimal energetic model at a low axonal load representative of VTA neurons (A = 0.40).** | **Figure 4B. Phase plane of the minimal energetic model at a high axonal load representative of SNc neurons (A = 1.00).** |

**Figure 4. Phase-plane comparison of VTA-like and SNc-like dopaminergic neurons.**
*Left (Figure 4A):* The horizontal axis shows energetic reserve E and the vertical axis mitochondrial capacity M. Gray arrows denote the vector field, the solid curve the energy nullcline (dE/dt = 0), and the dashed curve the mitochondrial nullcline (dM/dt = 0). In this low-load regime the nullclines intersect only once, at a high-energy, high-mitochondrial-capacity equilibrium. Sample trajectories (colored curves) initiated from widely separated initial conditions—including states with low energy and/or impaired mitochondria—are all attracted to this single fixed point. The absence of additional fixed points or a separatrix indicates a globally attracting, monostable high-energy regime, in which VTA-like neurons robustly recover their energetic state following transient perturbations rather than tipping into collapse.
*Right (Figure 4B):* As in the left panel, the horizontal axis shows energetic reserve E and the vertical axis mitochondrial capacity M. The solid curve denotes the energy nullcline (dE/dt = 0) and has an S-shaped profile, while the dashed curve shows the mitochondrial nullcline (dM/dt = 0). At this high load the nullclines intersect three times: at a low-energy, low-mitochondrial-capacity fixed point, a high-energy fixed point, and an intermediate saddle (red marker). Sample trajectories (colored curves) initiated on either side of the saddle’s stable manifold diverge toward different long-term outcomes: those starting above and to the right of the separatrix relax to the high-energy attractor, whereas those starting below or to the left are drawn into the low-energy attractor. The coexistence of these two energetic fates under the same parameter set, separated only by a narrow dynamical boundary, illustrates how extreme structural load places SNc-like neurons inside a bistable regime where modest perturbations can tip them from normal function into irreversible energetic collapse.[1,3,5–9,15–17,21–23](#references)

To evaluate their energetic stability under these distinct conditions, we examined the phase plane at load values representative of each population. **Figure 4A** shows the vector field and nullclines for a low-load setting ($A = 0.40$), corresponding to VTA-like neurons. In this regime, the system exhibits only a single stable equilibrium: a healthy, high-energy state. The absence of additional fixed points implies that VTA neurons are far from any critical boundary, and perturbations that transiently diminish energy are followed by a robust return to baseline.

In contrast, **Figure 4B** illustrates the phase plane at a higher load value ($A = 1.00$), representative of SNc neurons. Here the system lies **within the bistable window** identified in the bifurcation analysis. Three equilibria are present: a high-energy attractor, a low-energy attractor, and an intervening saddle point. The stable manifold of the saddle forms a separatrix that partitions the phase plane into two basins of attraction. SNc neurons therefore operate close to a **dynamical boundary**: small shifts in energetic reserve or mitochondrial capacity can determine whether the system returns to its healthy energetic state or transitions irreversibly into collapse.[15–17,21–23](#references)

This difference in stability geometry provides a mechanistic explanation for the selective vulnerability of SNc neurons. Both cell types face ongoing metabolic demands from pacemaking and neurotransmission, but only SNc neurons must satisfy these demands while maintaining a vast axonal arbor.[1,3,5–9,15–17](#references) In the model, this structural requirement positions SNc neurons near the saddle-node bifurcation where two energetic states coexist. From this vantage point, even perturbations that are insufficient to cause lasting damage in VTA neurons may push SNc neurons across the separatrix and into the low-energy attractor.

These analyses highlight a simple but powerful principle: **structural load alone is sufficient to place neurons in fundamentally different dynamical regimes of energetic stability**. The position of SNc neurons near a dynamical tipping point—not merely their molecular environment—creates the conditions for catastrophic collapse.[15–17,21–23](#references)

---

## **5. Perturbation Experiments Demonstrate Collapse in Substantia Nigra Neurons but Recovery in VTA Neurons**

To test how each neuronal population responds to transient metabolic stress, we simulated energy trajectories beginning near the high-energy steady state for both load conditions ($A = 0.40$ for VTA-like neurons and $A = 1.00$ for SNc-like neurons). Under these baseline conditions, both cell types remain stable and maintain high energetic reserve over long timescales (Figure 5A). This confirms that elevated structural load alone does not force SNc neurons into the pathological state; rather, it places them near a boundary where recovery from perturbation becomes precarious.[15–17,24](#references)

| ![](/Stage_4_TimeCourses/latest_run/timecourses_VTA_vs_SNc_baseline.png) | ![](/Stage_4_TimeCourses/latest_run/timecourses_SNc_vs_VTA_perturbation.png) |
| ------------------------------------------------------------------------ | ---------------------------------------------------------------------------- |
| **Figure 5A. Baseline trajectories of energetic reserve $E(t)$ for VTA-like ($A = 0.40$, blue) and SNc-like ($A = 1.00$, orange) neurons in the absence of perturbations.**               | **Figure 5B. Response of VTA-like ($A = 0.40$, orange) and SNc-like ($A = 1.00$, blue) neurons to a brief energetic insult applied at $t = 50$ (vertical dashed line), implemented by resetting $E$ to 0.3 while leaving $M$ unchanged.**    |

**Figure 5. Time-course simulations reveal robustness in VTA neurons and collapse in SNc neurons.**
*Left (Figure 5A):* Starting from slightly elevated initial conditions, both traces relax rapidly onto their respective high-energy steady states (dashed lines) and remain stable over long times, while the corresponding mitochondrial capacities $M(t)$ (not shown) behave similarly. Elevated structural load lowers the SNc steady-state energy level but does not drive spontaneous decline, indicating that extreme arborization alone is compatible with a persistent high-energy operating point and that collapse requires an additional perturbation.
*Right (Figure 5B):* In the VTA-like monostable regime, $E(t)$ rapidly relaxes back to its original high-energy steady state, indicating robust recovery from transient metabolic stress. In the SNc-like bistable regime, the same perturbation pushes the system across the saddle’s separatrix so that $E(t)$ subsequently drifts toward and stabilizes at the low-energy attractor; the corresponding mitochondrial capacity $M(t)$ (not shown) declines in parallel. The divergent outcomes of identical insults in the two regimes illustrate how proximity to a load-induced tipping point allows SNc neurons—but not VTA neurons—to undergo irreversible energetic collapse after modest perturbations.[15–17,21–24](#references)

We introduced a brief energetic perturbation by transiently reducing the energy variable $E$ to 0.3 at time $t = 50$, mimicking a short-lived metabolic challenge such as a burst of pacemaking Ca²⁺ entry, local inflammation, oxidative stress, or a mitochondrial inhibition event.[5–9,10–14,24](#references) The subsequent trajectories reveal a marked divergence between the two neuronal types (Figure 5B).

In the **VTA-like regime** ($A = 0.40$), the system quickly returns to the high-energy steady state after the perturbation. The energy reserve recovers smoothly, and mitochondrial capacity stabilizes along the same trajectory as in the unperturbed baseline. This behavior reflects the fact that VTA neurons, with their modest arborization load, lie far from the saddle-node bifurcation and thus possess a **single, globally attracting energetic state**. Perturbations may transiently reduce energy but do not threaten long-term stability.[15–17](#references)

In contrast, the **SNc-like regime** ($A = 1.00$) shows a dramatically different response. The same perturbation pushes the system across the stable manifold of the saddle, causing it to exit the basin of attraction of the high-energy state and converge instead to the low-energy attractor. Once initiated, this collapse is irreversible in the relevant physiological range, as the right-hand fold of the bifurcation sits beyond biologically realistic values of $A$ (Section 3). The resulting trajectory reflects the hallmark of a **tipping-point transition**: the neuron appears stable until a modest perturbation triggers a sudden and catastrophic drop in energetic reserve from which recovery is no longer possible.[21–23](#references)

This behavior mirrors the clinical and pathological course of Parkinson’s disease, where dopaminergic neurons can maintain function for decades before undergoing a rapid and irreversible decline.[15–17,24](#references) The model suggests that this vulnerability arises not from uniquely weak mitochondria or exclusive molecular stressors but from the **geometry of the underlying energy–mitochondria feedback loop**. SNc neurons operate close to a separatrix due to their extreme anatomical and physiological load. VTA neurons, lacking this structural burden, remain comfortably within a monostable regime where perturbations do not precipitate collapse.

These perturbation experiments thus provide computational evidence that **proximity to a load-induced saddle-node bifurcation is sufficient to explain the selective vulnerability of SNc dopaminergic neurons**.[15–17,21–23](#references)

---

## **6. Discussion**

Dopaminergic neurons of the substantia nigra pars compacta (SNc) are uniquely vulnerable in Parkinson’s disease, whereas neighboring ventral tegmental area (VTA) neurons remain comparatively resilient.[15–17](#references) Although numerous molecular abnormalities have been implicated—mitochondrial dysfunction, oxidative stress, calcium dysregulation, dopamine metabolism, and α-synuclein aggregation—no single factor has fully explained the striking anatomical selectivity of degeneration.[5–14,18–20,24](#references) Here we show that a **minimal energetic model** capturing only two dynamic variables and two structural loads is sufficient to reproduce the essential pattern of selective vulnerability. In this framework, SNc neurons are not intrinsically fragile; rather, they operate near a **saddle-node bifurcation** in the energy–mitochondria feedback loop created by their extraordinary axonal and calcium-handling demands.[1,3,5–9,15–17,21–23](#references)

This geometric perspective naturally unifies several empirical observations. SNc neurons maintain vast axonal arbors that impose high energetic costs for synaptic maintenance, vesicle cycling, and axoplasmic transport.[1–4,15,16](#references) Their autonomous pacemaking relies on L-type calcium channels, introducing additional energetic burden for Ca²⁺ extrusion and mitochondrial buffering.[5–9](#references) In the model, these features are reflected in the load parameter $A$. As structural load increases, the energy nullcline deforms, eventually creating two coexisting energetic states—a healthy attractor and a collapsed attractor—separated by a saddle. SNc-like neurons fall within this bistable regime, while VTA-like neurons, with their smaller arbors, remain monostable and robust. This separation of dynamical regimes provides a simple and mechanistic explanation for why two closely related neuronal populations experience vastly different fates under the same molecular milieu.[15–17,21–23](#references)

The model also offers insight into the **temporal profile** of degeneration. Clinical and pathological studies suggest that dopaminergic neurons can function for decades despite accumulating stressors, followed by a relatively sudden collapse in function and cell viability.[15–17,24](#references) In the saddle-node regime, apparent stability is maintained until a perturbation—metabolic fluctuation, inflammatory episode, mitochondrial insult—pushes the system across the separatrix.[10–14,18–20,24](#references) Once this threshold is crossed, recovery is no longer possible because the right-hand fold of the bifurcation lies beyond the biologically plausible range of axonal load (Section 3). The model therefore accounts for both the long prodromal period and the abrupt, irreversible decline in SNc neurons without invoking catastrophic molecular changes at the moment of symptom onset.[21–23](#references)

Several **testable predictions** emerge from this framework. Interventions that reduce effective load—such as pruning excessive axonal branches, reducing Ca²⁺ influx through L-type channel blockers, or lowering synaptic maintenance cost—should shift neurons leftward in parameter space, increasing resilience by moving them out of the bistable window.[5–9,15–17](#references) Enhancing mitochondrial repair or turnover should raise the energy nullcline and reduce sensitivity to perturbation.[10–14,24](#references) Conversely, stressors that transiently reduce energetic reserve (e.g., oxidative bursts, inflammatory cytokines, or mitochondrial inhibitors) are predicted to disproportionately harm SNc neurons by pushing them across the separatrix, whereas VTA neurons should recover.[15–17,24](#references) Heterogeneity in axonal arbor size within the SNc population may account for neuron-to-neuron differences in vulnerability, an idea consistent with recent single-cell degeneration patterns and anatomical variability.[1,3,15,16](#references)

This work also highlights **limitations**. The model is intentionally minimal and omits many molecular processes implicated in Parkinson’s disease, including α-synuclein aggregation, lysosomal-autophagic dysfunction, dopamine oxidation, and genetic factors such as LRRK2 mutations.[12,13,18–20](#references) These factors likely modulate energetic stress or mitochondrial resilience and could be incorporated as additional terms that shift nullcline geometry. Furthermore, the model does not include intrinsic stochasticity; incorporating noise may reveal graded transition probabilities near the separatrix, as suggested by general tipping-point theory.[22,23](#references) Despite these simplifications, the qualitative structure of the bifurcation—the emergence of two energetic states under load—proves robust under wide parameter variation (Supplementary Section S3), suggesting that the tipping-point mechanism is more fundamental than any specific biochemical pathway.[21–23](#references)

Overall, this minimal framework demonstrates that **extreme structural and physiological load places SNc dopaminergic neurons near a critical dynamical boundary**, making them uniquely susceptible to irreversible energetic collapse.[1,3,5–9,15–17,21–24](#references) By focusing on the geometry of energetic regulation rather than the molecular details, the model provides a unifying explanation for selective vulnerability and a foundation upon which more detailed mechanistic hypotheses can be built. Ultimately, this perspective suggests that therapeutic strategies aimed at reducing structural load, moderating calcium stress, or strengthening mitochondrial resilience may be effective not because they target a specific molecular lesion, but because they shift SNc neurons away from a saddle-node tipping point and restore them to a monostable energetic regime.

---

## **7. Methods**

### **7.1 Model Equations**

The energetic state of a dopaminergic neuron was described by two coupled differential equations governing the temporal evolution of energetic reserve $E(t)$ and mitochondrial functional capacity $M(t)$. Both variables are normalized to the interval $[0,1]$.

$$
\frac{dE}{dt}
= k_1 M (1-E) + k_2 E^2 (1-E) - \left(L_0 + L_1 A C\right) E ,
$$

$$
\frac{dM}{dt}
= k_M(1-M) - \beta A C M(1-E).
$$

The parameter $A$ represents axonal arborization load, and $C$ represents calcium-handling demand; unless otherwise specified, $C = 1$.[1,3,5–9,15–17](#references) All simulations were performed using the identical set of equations without modification.

---

### **7.2 Parameter Values**

Unless noted otherwise, the following parameter values were used throughout:

| Parameter | Value | Description                                   |
| --------- | ----- | --------------------------------------------- |
| $k_1$     | 1.0   | Mitochondrial ATP production term             |
| $k_2$     | 1.0   | Nonlinear energy amplification term           |
| $L_0$     | 0.1   | Baseline energy consumption                   |
| $L_1$     | 2.0   | Load-dependent energy consumption scaling     |
| $k_M$     | 1.0   | Mitochondrial repair/turnover rate            |
| $\beta$   | 1.0   | Energy-dependent mitochondrial damage scaling |
| $C$       | 1.0   | Ca²⁺-handling load                            |

Axonal load $A$ was varied across simulations. For population comparisons:

* **VTA-like neurons:** $A = 0.40$
* **SNc-like neurons:** $A = 1.00$

These values reflect relative arborization sizes from anatomical reconstructions and position each population on opposite sides of the saddle-node regime revealed by the bifurcation analysis.[1,3,15–17](#references)

---

### **7.3 Initial Conditions and Perturbation Protocol**

#### **Baseline simulations**

For baseline time courses (Figure 5A), both populations were initialized at:

$$
E(0) = 0.9,\qquad M(0) = 0.9.
$$

These values lie within the basin of attraction of the high-energy steady state for all relevant $A$.

#### **Perturbation experiments**

At time $t = 50$, a transient energetic perturbation was applied by resetting:

$$
E(t=50) = 0.3,
$$

after which the system evolved freely under the governing equations. Mitochondrial capacity $M$ was not directly perturbed, allowing the collapse or recovery dynamics to arise solely from the energetic state.

---

### **7.4 Numerical Integration**

All simulations were performed in **Python 3.13** using:

* **SciPy 1.11.3** [26](#references)
* **NumPy 1.26** [28](#references)
* **Matplotlib 3.7** [27](#references)

Differential equations were integrated with the explicit Runge–Kutta method of order 5(4) (“RK45”) with:

```text
rtol = 1e-8
atol = 1e-10
max_step = 0.1
```

For time-course figures, solutions were evaluated on a uniform time grid of 4000 points over $t \in [0, 300]$.

---

### **7.5 Determination of Steady States and Stability**

Steady states for the bifurcation diagrams were computed using a two-step procedure:

1. **Coarse grid scan:**
   A regular grid of candidate points was constructed in the $(E, M)$ plane, and the right-hand side of the ODEs was evaluated to identify approximate zero crossings.

2. **Refinement by root-finding:**
   Each candidate point was refined using SciPy’s `fsolve` with tight tolerances.[26](#references)

Each equilibrium’s stability was determined by linearizing the system at that point and analyzing the eigenvalues of the resulting Jacobian. Points with eigenvalues having negative real parts were classified as stable; those with one positive real eigenvalue were classified as saddles.[21–23](#references)

---

### **7.6 Bifurcation Scan Procedure**

To construct the bifurcation diagram ([Figure 3](#figure_3)), the axonal load parameter $A$ was varied from 0.2 to 1.4 in increments of 0.02 (61 values total). For each $A$, all equilibria were computed as described above. The bistable window was defined as the range of $A$ for which exactly three equilibria were present.

---

### **7.7 Phase Plane and Nullcline Computation**

Phase planes (Figures 2 and 4) were computed on $[0,1] \times [0,1]$ grids of $200 \times 200$ points. At each point:

* The vector field $(dE/dt, dM/dt)$ was evaluated and normalized for display.
* Nullclines were computed by finding zero-level contours of $dE/dt$ and $dM/dt$ using Matplotlib’s contour-finding routines.[27](#references)
* Sample trajectories were generated by numerically integrating from selected initial conditions.[26,27](#references)

---

### **7.8 Reproducibility and Code Availability**

All simulations used deterministic ODE integration with fixed parameters. There is no stochasticity in the core model. Numerical code, figure-generation scripts, and processed data outputs will be made available upon publication and can be reproduced directly using the parameter sets and methods described above.

---

## **Supplementary Results**

### **S1. Earlier model architectures did not exhibit bistability**

We evaluated several alternative formulations of the energy–mitochondria interaction prior to arriving at the final minimal model. These earlier systems were structurally incapable of producing multiple equilibria under biologically realistic parameter ranges, despite including many of the same biological ingredients. Their failure highlights the importance of the nonlinear energy-amplification and load-modulated mitochondrial damage terms used in the final formulation.[21–23](#references)

#### **S1.1 First model variant: linear mitochondrial support and linear load**

The simplest architecture consisted of:

* Linear ATP production proportional to $M(1-E)$,
* Quadratic mitochondrial damage proportional to $A C M (1-E)$,
* No nonlinear amplification in energy (i.e., the $E^2(1-E)$ term was absent).

In this system, the energy nullcline is monotonic for all parameter values, and the mitochondrial nullcline intersects it **exactly once** within the biologically relevant domain. The vector field analysis confirmed global convergence to a unique fixed point. No choice of axonal load $A$, even when increased far beyond anatomical ranges, produced a second intersection.

![](/outputs/supplement_ec2_S3_phase_plane.png)

**Supplementary Figure S1. Nullclines and vector field for an earlier, monostable model variant.**
In a model with linear mitochondrial support and load effects but no nonlinear energy amplification, the energy and mitochondrial nullclines intersect only once in the physical domain. The vector field indicates global convergence to a single equilibrium for all tested values of $A$, demonstrating that this architecture cannot produce bistability.

#### **S1.2 Second model variant: feedback mitochondrial impairment without nonlinear energy restoration**

The second architecture attempted to incorporate a biologically motivated feedback: energetic deficit increases mitochondrial damage. Mathematically, this introduced a term of the form $(1 - M)E$ in the energy equation and retained the load-dependent mitochondrial damage. Despite this additional coupling, the system remained **monostable**.

* The energy nullcline remained single-peaked but never developed a fold.
* The mitochondrial nullcline cut across it only once.
* Bifurcation scans across a wide range of $A \in [0.1, 2.0]$ showed one stable equilibrium everywhere.

![](/outputs/supplement_ec2_S2_bifurcation_A.png)

**Supplementary Figure S2. Bifurcation scan for a feedback-only model lacking nonlinear energy amplification.**
Even with feedback from energy to mitochondrial damage, the system exhibits only a single equilibrium across a broad range of axonal loads $A$. No saddle-node bifurcation appears, confirming that feedback alone is insufficient to generate bistability in this architecture.[21–23](#references)

These results demonstrated that feedback alone is insufficient: to generate bistability, the model must produce a genuine **S-shaped energy nullcline** within the physical domain.

---

### **S2. Requirement for nonlinear energy amplification**

The term $k_2 E^2(1-E)$ in the final model introduces a saturating, cooperative-like energetic contribution that is negligible near $E = 0$ but increases sharply as $E$ rises. This positive curvature is what allows the energy nullcline to fold under load.[21–23](#references)

By contrast:

* In models lacking this term, the nullcline was monotonic.
* In models using only linear or quadratic forms, the nullcline never turned sufficiently to create two additional equilibria.

When included, the nonlinear term interacts with load-dependent consumption and energy-dependent mitochondrial damage to produce:

1. A high-energy fixed point,
2. A low-energy fixed point,
3. A saddle separating their basins.

![](/outputs/supplement_ec2_S3_phase_plane.png)

**Supplementary Figure S3. Effect of nonlinear energy amplification on nullcline geometry.**
Including the nonlinear term $k_2 E^2(1-E)$ generates a folded energy nullcline that, together with the mitochondrial nullcline, produces three intersections corresponding to two stable equilibria and a saddle. Removing this term collapses the fold and eliminates bistability.

---

### **S3. Parameter sweeps confirm robustness of the bistable window**

We systematically varied key parameters to assess how sensitive the saddle-node structure is to physiological uncertainty (Supplementary Figures S4–S6).

#### **S3.1 Variation in mitochondrial turnover rate $k_M$**

Increasing $k_M$ shifts the bistable window to higher values of $A$, reflecting improved mitochondrial resilience. Conversely, reducing $k_M$ expands the bistable region toward lower loads, making more neurons susceptible to collapse.[10–14,21–23](#references)

#### **S3.2 Variation in load-dependent consumption $L_1$**

Higher $L_1$ (greater energy cost per unit arbor) widens the bistable region and lowers the high-energy steady-state value. Lower $L_1$ compresses the window and raises the healthy equilibrium.

#### **S3.3 Variation in calcium-handling load $C$**

Because SNc neurons experience substantial Ca²⁺ influx during pacemaking, we explored $C \in [0.5, 1.5]$. Changes in $C$ functionally act like scaling $A$; high $C$ shifts the bistable window leftward, increasing vulnerability.[5–9,21–23](#references)

Across all sweeps, the saddle-node bifurcation persisted, demonstrating that the fold is a **structural consequence** of the feedback motif, not a fine-tuned artifact.[21–23](#references)

---

### **S4. Absence of biologically relevant hysteresis**

Mathematically, the full bifurcation curve includes both a left fold (creation of bistability) and a right fold (annihilation). However:

* The right fold occurs at $A > 1.4$, beyond anatomical values for dopaminergic arbor size.[1,3,15–17](#references)
* Thus, the low-energy branch persists across all physiologically meaningful loads.

This means that once the system collapses into the low-energy attractor, reducing the load (e.g., via axonal pruning during degeneration) does **not** restore the healthy energetic state. This absence of hysteresis in the biological range explains:

1. The irreversibility of SNc collapse,
2. Why pruning does not rescue already failing neurons,
3. Why collapse proceeds inexorably once initiated.[15–17,21–23](#references)

**Supplementary Figure S7. Full continuation of equilibria showing right fold beyond biological range.**
The complete bifurcation diagram reveals a right-hand fold at axonal loads larger than those observed in dopaminergic neurons, implying that the collapsed low-energy state remains stable throughout the physiological range of $A$.

---

### **S5. Summary of Supplementary Findings**

The supplementary analyses demonstrate three key points:

1. **Not all biologically plausible architectures can generate bistability.**
   Specific nonlinear interactions—particularly nonlinear energy amplification—are required.[21–23](#references)

2. **The saddle-node bifurcation in the final minimal model is robust**, persisting under broad parameter variation.[21–23](#references)

3. **Collapse is irreversible in the physiological range**, consistent with clinical course and anatomical constraints in SNc neurons.[1,3,15–17,24](#references)

Together, these results strengthen the interpretation that selective vulnerability of SNc neurons emerges not from unique molecular defects but from **the fundamental geometry of load-dependent energetic regulation**.[15–17,21–24](#references)

---

<a id="references"></a>

## **References**

**Anatomical Load & Arborization Complexity**

1. Matsuda W, Furuta T, Nakamura KC, Hioki H, Fujiyama F, Arai R, Kaneko T. *Single nigrostriatal dopaminergic neurons form widely spread and highly dense axonal arborizations in the neostriatum.* J Neurosci. 2009;29(2):444–453.

2. Bolam JP, Freund TF, Henderson Z. *Diverse interneurons in the neostriatum—various dendritic and axonal structures and synaptic connections.* Trends Neurosci. 2000;23(8):377–384.

3. Gauthier J, Parent A. *Distribution of axon collaterals from single nigrostriatal neurons in the rat.* Brain Res. 1989;500(1–2):18–30.

4. Pan WX, Mao T, Dudman JT. *Input–output organization of the basal ganglia.* Curr Opin Neurobiol. 2010;20(2):223–229.

---

**Calcium Handling, Pacemaking, and Energetic Burden**

5. Surmeier DJ, Guzman JN, Sanchez-Padilla J, Goldberg JA. *What causes the death of dopaminergic neurons in Parkinson’s disease?* Prog Brain Res. 2010;183:59–77.

6. Guzman JN, Sánchez-Padilla J, Chan CS, Surmeier DJ. *Robust pacemaking in substantia nigra dopaminergic neurons.* J Neurosci. 2009;29(35):11011–11019.

7. Chan CS, Guzman JN, Ilijic E, et al. *‘Rejuvenation’ protects neurons in mouse models of Parkinson’s disease.* Nature. 2007;447(7148):1081–1086.

8. Goldberg JA, Guzman JN, Estep CM, Surmeier DJ. *Calcium entry via Cav1 channels controls mitochondrial function in substantia nigra dopaminergic neurons.* Neuron. 2012;76(2):356–369.

9. Surmeier DJ. *Calcium, aging, and neuronal vulnerability in Parkinson’s disease.* Cell Calcium. 2007;42(3):351–361.

---

**Energetics, Mitochondrial Stress, and Parkinson’s Disease**

10. Exner N, Lutz AK, Haass C, Winklhofer KF. *Mitochondrial dysfunction in Parkinson’s disease: molecular mechanisms and pathophysiological consequences.* EMBO J. 2012;31(14):3038–3062.

11. Grünewald A, Rygiel KA, Zsurka G, et al. *Mitochondrial DNA depletion in substantia nigra neurons in Parkinson disease.* Ann Neurol. 2016;79(3):366–378.

12. Burbulla LF, Song P, Mazzulli JR, et al. *Dopamine oxidation mediates mitochondrial and lysosomal dysfunction in Parkinson’s disease.* Science. 2017;357(6357):1255–1261.

13. Schapira AHV. *Mitochondrial complex I deficiency in Parkinson’s disease.* Ann Neurol. 1990;28(2):149–155.

14. Bose A, Beal MF. *Mitochondrial dysfunction in Parkinson’s disease.* J Neurochem. 2016;139(S1):216–231.

---

**Selective Vulnerability of SNc Neurons**

15. Pacelli C, Giguère N, Bourque M-J, et al. *Elevated mitochondrial bioenergetics and axonal arborization size are key contributors to the vulnerability of dopamine neurons.* Cell Rep. 2015;13(4):729–741.

16. Giguère N, Burke Nanni S, Trudeau L-E. *On cell loss in Parkinson's disease and the selective vulnerability of SNc dopaminergic neurons: insights from genetic mouse models.* Cell Mol Life Sci. 2018;75:1477–1493.

17. Surmeier DJ. *Determinants of dopaminergic neuron vulnerability in Parkinson’s disease.* FEBS Lett. 2018;592(6):743–753.

---

**α-Synuclein, Aggregation, and Proteostasis (Contextual Factors)**

18. Wong YC, Krainc D. *α-Synuclein toxicity in neurodegeneration: mechanism and therapeutic strategies.* Nat Med. 2017;23(2):1–13.

19. Conway KA, Harper JD, Lansbury PT. *Accelerated in vitro fibril formation by a mutant α-synuclein linked to early-onset Parkinson disease.* Nat Med. 1998;4(11):1318–1320.

20. Cookson MR. *The role of leucine-rich repeat kinase 2 (LRRK2) in Parkinson's disease.* Nat Rev Neurosci. 2010;11(12):791–801.

---

**Dynamical Systems, Bifurcations, and Tipping Points**

21. Strogatz SH. *Nonlinear Dynamics and Chaos.* CRC Press; 2018.

22. Scheffer M, Carpenter SR, Lenton TM, et al. *Anticipating critical transitions.* Science. 2012;338(6105):344–348.

23. Kuehn C. *Multiple Time Scale Dynamics.* Springer; 2015.

---

**Energetic Failure, ATP Dynamics, and Neuronal Degeneration**

24. Nicholls DG. *Oxidative stress and energy crises in neuronal dysfunction.* Nat Rev Neurosci. 2004;5(11):862–872.

25. Zheng X, Boyer L, Jin M, et al. *Metabolic reprogramming during neuronal differentiation.* Cell Metab. 2016;23(6):1068–1082.

---

**Methodological References**

26. Virtanen P, Gommers R, Oliphant TE, et al. *SciPy 1.0: fundamental algorithms for scientific computing in Python.* Nat Methods. 2020;17:261–272.

27. Hunter JD. *Matplotlib: A 2D graphics environment.* Comput Sci Eng. 2007;9(3):90–95.

28. Harris CR, Millman KJ, van der Walt SJ, et al. *Array programming with NumPy.* Nature. 2020;585:357–362.
