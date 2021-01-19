;limitations/Assumptions: NOTE: This is far from a comprehensive list. For a more-detailed list refer to our Overview, Design concepts, Details model description available at: and on the Lanzaslab github page (https://github.com/lanzaslab/droplet-ABM)
;; 1.) Movement in model is extremely simplified - In actuallity, infectious and healthy people will likely enter rooms through the same door, move through pre-designated pathways, and may interact more frequently than we represent here. Which provides opportunities for localized droplet or airborne droplet-expulsion in these locations.
;; 2.) No people transition from infected to infectious in this model.
;; 3.) People are not added to/removed from cohorts during droplet-expulsion simulations, but cohorts may be completely replaced.
;; 4.) People have free reign of the space in world (i.e., there's no obstruction by seating, tables, etc.).
;; 5.) Infection is simplified. The mechanism of infection is assumed to be the same for droplets of all sizes, and is simualted as inhalation.
;; 6.) Airflow is two-dimensional.

; Note: wrapping has been turned off for this model (e.g., turtles cannot move from patch (0, min(pycor)) to patch (0, max(pycor)) in one step)

globals [

  avg.dist              ; Tracks the average interpersonal distance at each tick.
  avg.exposureTime      ; Tracks the average number of timesteps when people were exposed to >= 1 infectious virion at each tick.
  avg.PatchInfectiousness ; Tracks the average probability that exposure to patch virions will lead to infection.
  beta_asymp            ; Variable used as part the function for generating a log-normal random draw for pathogen spread from ASYMPTOMATIC individuals.
  beta_cough_droplets      ; Variable used as part the function for generating a log-normal random draw for droplet creation from COUGHING.
  beta_speak_droplets      ; Variable used as part the function for generating a log-normal random draw for droplet creation from SPEAKING.
  beta_symp             ; Variable used as part the function for generating a log-normal random draw for pathogen spread from SYMPTOMATIC individuals.
  cohort                ; Keep track of what cohort is present.
  cohort-endTime        ; Denotes the tick value when cohorts be replaced.
  firstInfectTime       ; Records the tick at which the first successful transmission event occurs.
  lastInfectTime        ; Records the tick at which the last susceptible individual was infected (only) relevant if ALL susceptibles were infected.
  M_asymp               ; Variable used as part the function for generating a log-normal random draw for pathogen spread from ASYMPTOMATIC individuals.
  M_cough_droplets         ; Variable used as part the function for generating a log-normal random draw for droplet creation from COUGHING.
  M_speak_droplets         ; Variable used as part the function for generating a log-normal random draw for droplet creation from SPEAKING.
  M_symp                ; Variable used as part the function for generating a log-normal random draw for pathogen spread from SYMPTOMATIC individuals.
  num_asymptomatic      ; The number of asymptomatic "infectious" people that spawn in the first cohort. The probability of an infectious agent being asymptomatic is (n_infectious *(1 - symp-pr)).
  num_symptomatic       ; The number of symptomatic "infectious" people that spawn in the first cohort. The probability of an infectious agent being symptomatic is (n_infectious * symp-pr).
  num_completelySusceptible ; Counts the number of individuals in a cohort that are completely susceptible to infection.
  num_reducedSusceptible  ; Counts the number of individuals in a cohort with reduced susceptibility to infection.
  patchContamination.list  ; List of the number of contaminated patches present at each time step.
  patchInfectiousness.list ; List of the mean patch infectiousness values observed throughout the simulation.
  personDist.list       ; List of average distance between people at each tick. Note: This is a secondary output of interest.
  S_asymp               ; Variable used as part the function for generating a log-normal random draw for pathogen spread from ASYMPTOMATIC individuals.
  S_cough_droplets         ; Variable used as part the function for generating a log-normal random draw for droplet creation from COUGHING.
  S_speak_droplets         ; Variable used as part the function for generating a log-normal random draw for droplet creation from SPEAKING.
  S_symp                ; Variable used as part the function for generating a log-normal random draw for pathogen spread from SYMPTOMATIC individuals.
  totalInfected         ; Tracks the total number of infected people over the course of the simulation (i.e., running sum). Note: This is the primary output of interest.
  totalInfected.list    ; List of the total number of infected people at each tick. Note: This is a primary output of interest.
  Vt_diam3              ; The terminal velocity (m/s) of a respiratory droplet with a 3-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam6              ; The terminal velocity (m/s) of a respiratory droplet with a 6-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam12             ; The terminal velocity (m/s) of a respiratory droplet with a 12-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam20             ; The terminal velocity (m/s) of a respiratory droplet with a 20-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam28             ; The terminal velocity (m/s) of a respiratory droplet with a 28-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam36             ; The terminal velocity (m/s) of a respiratory droplet with a 36-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam45             ; The terminal velocity (m/s) of a respiratory droplet with a 45-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam62.5           ; The terminal velocity (m/s) of a respiratory droplet with a 62.5-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam87.5           ; The terminal velocity (m/s) of a respiratory droplet with a 87.5-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam112.5          ; The terminal velocity (m/s) of a respiratory droplet with a 112.5-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam137.5          ; The terminal velocity (m/s) of a respiratory droplet with a 137.5-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam175            ; The terminal velocity (m/s) of a respiratory droplet with a 175-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam225            ; The terminal velocity (m/s) of a respiratory droplet with a 225-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam375            ; The terminal velocity (m/s) of a respiratory droplet with a 375-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  Vt_diam750            ; The terminal velocity (m/s) of a respiratory droplet with a 750-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).

  ; The following global variables are ones defined on the interface tab

  ;; choirSim           ; Logical. Describes if we want to simulate the circumstances of the choir prcatice SARS-CoV2 case study whereby 51 people were likely infected by one infectious individual at the practice (Hamner et al. 2020). The primary effect of this parameter is that cohort members are rearranged 3 times after different lengths of times (i.e., 40 ticks, 50 ticks, 15 ticks). This simulates, as best we can, the movement of individuals within an between rooms of the church where they were practicing. Note that in the actual case, some individuals separated to practice in a separate room for 50 mins. There is no information available about the number or specifics of people that moved, or the location of the separate room relative to the primary rehersal space. Additionally, the final rearrangement (i.e., after 40 + 50 + 15 ticks) puts indivivudals back into their original locations to simulate them moving back into their original seats.
  ;; cohort-dur         ; The number of ticks (minutes) that each cohort lasts (i.e., how long each people cohort spends in the simulation)
  ;; cough_airflow-angle ; Controls the angle of airflow associated with coughing events. Affects the spread of droplets during a given droplet-expulsion event originating from symptomatic individuals.
  ;; cough-frequency    ; Probability that cones of infection stemming from symptomatic individuals will be parameterized using cough_airflow-angle, cough_spread-dist.mean, and cough_spread-dist.sd, instead of the asymptomatic counterparts.
  ;; cough_spread-dist.mean   ; The mean distance from symptomatic infectious people that droplets may be spread when coughing. This value will be used to generate a lognormal distribution from which droplet-expulsion-spread distance will be randomly drawn when a symptomatic agent triggers a droplet-expulsion event.
  ;; cough_spread-dist.sd     ; The standard deviation distance, given a cough_spread-dist.mean value, from symptomatic infectious people that droplets may be spread when coughing. This value will be used to generate a lognormal distribution from which droplet-expulsion-spread distance will be randomly drawn
  ;; diffusionRate      ; The rate (m^2 / min) at which droplets spread to adjacent patches. For simplicity we assume a standardized rate for all droplet sizes. Additionally, note that we assume droplets move in all directions at the same rate (i.e., the unit here can also be considered to be m^3 / min). However, droplets may only only move in two dimensions (i.e., they neither rise above the the patch "ceiling" or settle on the "ground" to explicitly become inactivated. Instead, dropletDecay is the parameter that controls the proportion of droplets removed from patches.
  ;; dropletDecay       ; Droplet decay rate per minute.
  ;; expectorate-height ; The height (in m) at which droplets are expelled. This is also the maximum vertical height of the simulated world, and the height used to in area volume calculations. (i.e., patch volumes are 1 m X 1 m X expectorate-height m).
  ;; evalQuanta         ; Logical. If true, all the simulation does is distribute droplets (number informed by speaking mean and sd) with diameter <= 20 across all patches, run an infection check , and report the proportion of n that is infected. The goal is tobe used in the determination of what constitutes a "quantum" (i.e., a number of airborne droplets suficient to infect (1 - 1/e) * 100% of the population).
  ;; face-northward       ; Logical variable that controls whether people only look northward within a range of 90 degrees (if TRUE) or face a random direction (if FALSE).
  ;; grid-height        ; The number of rows present in the grid representing the room in which agents interact. Note: cells in the matrix (i.e., patches), regardless of how many there are, represent 1 m X 1 m areas.
  ;; grid-width         ; The number of columns present in the grid representing the room in which agents interact. Note: cells in the matrix (i.e., patches), regardless of how many there are, represent 1 m X 1 m areas.
  ;; maskRisk-mod       ; The proportion of droplets exhaled/inhaled by mask wearers. Ex. if maskRisk-mod is 0.1, infectious will only expectorate 10% of the infectious droplets they would otherwise expel, and susceptible individuals will only inhale 10% of the droplets they would normally inhale.
  ;; mod_group          ; Controls what agent variables are modified by risk mod. Takes the values "sus," "inf," or "sus_inf" (representing susceptible agents only, infectious agents only, or both). If "sus," only exposure-risk is adjusted. If "inf," only expectorate-risk is adjusted. If "sus_inf," both of these variables are updated.
  ;; mod-proportion     ; Describes the probability that susceptible individuals will have their infection probability modified by the maskRisk-mod parameter. This parameter is used to vary the proportion of individuals minimizing their disease risk in the population (e.g., through the use of PPE).
  ;; n                  ; The total number of (i.e., both "healthy" and "infectious") people turtles that spawn in each cohort.
  ;; n_infectious       ; The number of infectious people turtles that spawn in the first cohort. This is a subset of n, not additional turtles. Infectious people default to "asymptomatic" status.
  ;; num-cohorts        ; The number of people cohorts observed during the simulation. Note: all cohorts interact with the same grid (i.e., world), but do not exist in the world at the same time
  ;; numReturnVents     ; The number of patches designated as return vents. Cannot exceed grid-width value if ventilReturnWall is one of "north," "south," "up," or "down." Cannot exceed grid-height value if ventilReturnWall is one of "east," "west," "left," or "right."
  ;; numSupplyVents     ; The number of patches designated as supply vents. Cannot exceed grid-width value if ventilSupplyWall is one of "north," "south," "up," or "down." Cannot exceed grid-height value if ventilReturnWall is one of "east," "west," "left," or "right."
  ;; personPerPatch-cap ; The maximum number of people that may possibly exist within a single patch.
  ;; rearrange-cohort   ; Logical variable describing whether or not the cohort-replace effectively becomes a rough proxy for movement of people within the room. If TRUE, the first "cohort" never leaves the room, rather, they are re-distributed according to the social-distance and personPerPatch-cap values set.
  ;; showArrows         ; Logical. If TRUE, airArrows will be visible. If FALSE, they will be hidden.
  ;; social-distance    ; The interpersonal distance, in meters, that agents seek to maintain over the course of the simulation.
  ;; speak_airflow-angle ; The angle of airflow associated with breathing events. Affects the spread of droplets during a given droplet-expulsion event from asymptomatic individuals.
  ;; speak_spread-dist.mean   ; The mean distance from asymptomatic infectious people that droplets may be spread when breathing. This value will be used to generate a lognormal distribution from which droplet-expulsion-spread distance will be randomly drawn when an asymptomatic agent triggers a droplet-expulsion event.
  ;; speak_spread-dist.sd     ; The standard deviation distance, given an speak_spread-dist.mean value, from asymptomatic infectious people that droplets may be spread when breathing. This value will be used to generate a lognormal distribution from which droplet-expulsion-spread distance will be randomly drawn when an asymptomatic agent triggers a droplet-expulsion event.
  ;; symp-pr             ; Probability that infectious people will be "symptomatic," instead of having the default "asymptomatic" status.
  ;; ventilation         ; Logical variable describing if airflow will move droplets throughout patches during the simulation.
  ;; ventil_movementRate ; Describes the proportion of droplets suspended in the air that will move to another patch (directed by location of return vent(s)) at each tick (Note: this can also be considered to be the "percent indoor air replacement / minute"). Note that this value must be <= 1 (i.e., 100%).
  ;; ventil_removalRate ; Describes the proportion of droplets on return vent patch(es) that will be removed from the simulation due to filtration.
  ;; ventilReturnWall  ; Takes one value "north," "south," "east," "west," OR "up," "down," "right," "left." Describes the wall of the simulated world that return vents will be located on.
  ;; ventilSupplyWall  ; Takes one value "north," "south," "east," "west," OR "up," "down," "right," "left." Describes the wall of the simulated world that supply vents will be located on.
  ;; virionsPerML       ; Number of virions per mL of droplet fluid.
  ;; virionRisk         ; The risk of infection given exposure to a single virion.
  ;; vol_B            ; The rate of air inhaled by individuals in patches (in cubic meters per minute).
]

breed [people person]   ; infectious, exposed, and susceptible agents in this model are people
breed [airArrows airArrow]  ; airArrows are used to control the direction of airflow-induced droplet movement in the simulation.

patches-own [

  additionalDroplets_size3         ; Used in the move-air action to add a number of droplets in a size class with a mean size of 3 micrometers.
  additionalDroplets_size6         ; Used in the move-air action to add a number of droplets in a size class with a mean size of 6 micrometers.
  additionalDroplets_size12        ; Used in the move-air action to add a number of droplets in a size class with a mean size of 12 micrometers.
  additionalDroplets_size20        ; Used in the move-air action to add a number of droplets in a size class with a mean size of 20 micrometers.
  additionalDroplets_size28        ; Used in the move-air action to add a number of droplets in a size class with a mean size of 28 micrometers.
  additionalDroplets_size36        ; Used in the move-air action to add a number of droplets in a size class with a mean size of 36 micrometers.
  additionalDroplets_size45        ; Used in the move-air action to add a number of droplets in a size class with a mean size of 45 micrometers.
  additionalDroplets_size62.5      ; Used in the move-air action to add a number of droplets in a size class with a mean size of 62.5 micrometers.
  additionalDroplets_size87.5      ; Used in the move-air action to add a droplets in a size class with a mean size of 87.5 micrometers.
  additionalDroplets_size112.5     ; Used in the move-air action to add a number of droplets in a size class with a mean size of 112.5 micrometers.
  additionalDroplets_size137.5     ; Used in the move-air action to add a number of droplets in a size class with a mean size of 137.5 micrometers.
  additionalDroplets_size175       ; Used in the move-air action to add a number of droplets in a size class with a mean size of 175 micrometers.
  additionalDroplets_size225       ; Used in the move-air action to add a number of droplets in a size class with a mean size of 225 micrometers.
  additionalDroplets_size375       ; Used in the move-air action to add a number of droplets in a size class with a mean size of 375 micrometers.
  additionalDroplets_size750       ; Used in the move-air action to add a number of droplets in a size class with a mean size of 750 micrometers.
  can-sprout?            ; Logical variable describing if patches are far enough away from those with people in them that new people can sprout a while keeping the desired social-distance value.
  person-count           ; Counts the number of people in the cell.
  droplets_size3         ; Counts the number of droplets in a size class with a mean size of 3 micrometers.
  droplets_size6         ; Counts the number of droplets in a size class with a mean size of 6 micrometers.
  droplets_size12        ; Counts the number of droplets in a size class with a mean size of 12 micrometers.
  droplets_size20        ; Counts the number of droplets in a size class with a mean size of 20 micrometers.
  droplets_size28        ; Counts the number of droplets in a size class with a mean size of 28 micrometers.
  droplets_size36        ; Counts the number of droplets in a size class with a mean size of 36 micrometers.
  droplets_size45        ; Counts the number of droplets in a size class with a mean size of 45 micrometers.
  droplets_size62.5      ; Counts the number of droplets in a size class with a mean size of 62.5 micrometers.
  droplets_size87.5      ; Counts the number of droplets in a size class with a mean size of 87.5 micrometers.
  droplets_size112.5     ; Counts the number of droplets in a size class with a mean size of 112.5 micrometers.
  droplets_size137.5     ; Counts the number of droplets in a size class with a mean size of 137.5 micrometers.
  droplets_size175       ; Counts the number of droplets in a size class with a mean size of 175 micrometers.
  droplets_size225       ; Counts the number of droplets in a size class with a mean size of 225 micrometers.
  droplets_size375       ; Counts the number of droplets in a size class with a mean size of 375 micrometers.
  droplets_size750       ; Counts the number of droplets in a size class with a mean size of 750 micrometers.
  returnVent             ; Denotes if patch is a return vent.
  supplyVent             ; Denotes if patch is a supply vent.
  totalDroplets          ; Counts the total number of droplets present in the patch.
  transmissionRisk       ; Tracks the probability that susceptible people on the patch will be infected on a given time point. This is the product of virionCount and virionRisk.
  virionCount            ; Counts the number of virions in the patch on a given time step.
]

people-own [

  droplets_size3AtInf         ; Counts the number of droplets in a size class with a mean size of 3 micrometers that were present in the containing patch when the individual was infected.
  droplets_size6AtInf         ; Counts the number of droplets in a size class with a mean size of 6 micrometers that were present in the containing patch when the individual was infected.
  droplets_size12AtInf        ; Counts the number of droplets in a size class with a mean size of 12 micrometers that were present in the containing patch when the individual was infected.
  droplets_size20AtInf        ; Counts the number of droplets in a size class with a mean size of 20 micrometers that were present in the containing patch when the individual was infected.
  droplets_size28AtInf        ; Counts the number of droplets in a size class with a mean size of 28 micrometers that were present in the containing patch when the individual was infected.
  droplets_size36AtInf        ; Counts the number of droplets in a size class with a mean size of 36 micrometers that were present in the containing patch when the individual was infected.
  droplets_size45AtInf        ; Counts the number of droplets in a size class with a mean size of 45 micrometers that were present in the containing patch when the individual was infected.
  droplets_size62.5AtInf      ; Counts the number of droplets in a size class with a mean size of 62.5 micrometers that were present in the containing patch when the individual was infected.
  droplets_size87.5AtInf      ; Counts the number of droplets in a size class with a mean size of 87.5 micrometers that were present in the containing patch when the individual was infected.
  droplets_size112.5AtInf     ; Counts the number of droplets in a size class with a mean size of 112.5 micrometers that were present in the containing patch when the individual was infected.
  droplets_size137.5AtInf     ; Counts the number of droplets in a size class with a mean size of 137.5 micrometers that were present in the containing patch when the individual was infected.
  droplets_size175AtInf       ; Counts the number of droplets in a size class with a mean size of 175 micrometers that were present in the containing patch when the individual was infected.
  droplets_size225AtInf       ; Counts the number of droplets in a size class with a mean size of 225 micrometers that were present in the containing patch when the individual was infected.
  droplets_size375AtInf       ; Counts the number of droplets in a size class with a mean size of 375 micrometers that were present in the containing patch when the individual was infected.
  droplets_size750AtInf       ; Counts the number of droplets in a size class with a mean size of 750 micrometers that were present in the containing patch when the individual was infected.
  cohort-person         ; Denotes what cohort the person belongs to.
  expectorate-risk      ; Denotes proportion of droplets agents expel on any given timestep. Defaults to 1. Changes if people are wearing masks.
  exposure-risk         ; Denotes agents' probability of infection given exposure to infectious agents. Defaults to 1 (i.e., complete susceptibility).
  infected?             ; Logical variable describing if people have been infected by contaminated patches.
  infectious?           ; Logical variable describing if people can spread the pathogen.
  mask?                 ; Logical variable describing if people are wearing a mask.
  rememberPatch         ; Sets a patch that people will remember, so that they may return to it later on.
  symptomatic?          ; Logical variable describing if people are coughing to spread the contagion. If TRUE, spread will be dictated by cough_airflow-angle, cough_spread-dist.mean, and cough_spread-dist.sd parameter values. If FALSE, but infectius? is TRUE, spread will be dictated by speak_airflow-angle, speak_spread-dist.mean, and speak_spread-dist.sd parameter values.
  timesteps_exposed     ; Running count of the number of timesteps when an agent was exposed to >= 1 infectious virion.

]

to setup

  if(n > (grid-height * grid-width * personPerPatch-cap))[ ; If there are too many individuals to fit in the simulated world.

    show "Not enough patches to hold n individuals given the current personPerPatch-cap. Simulation aborted." ; tell users why the simulation will not move forward.
    stop ; abort the simulation.

  ]

  if(ventilation = TRUE)[ ; if airflow is TRUE there are a few other things could cause errors. We prevent them here.

  if(ventilReturnWall = "north" OR ventilReturnWall = "south" OR ventilReturnWall = "up" OR ventilReturnWall = "down")[ ; determine whether or not there are enough patches to support vents

    if(numReturnVents > grid-width)[ ; If there are too many vents.

      show "numReturnVents exceeds the number of patches on the designated wall. Simulation aborted." ; tell users why the simulation will not move forward.
      stop ; abort the simulation.

    ]
  ]

  if(ventilSupplyWall = "north" OR ventilSupplyWall = "south" OR ventilSupplyWall = "up" OR ventilSupplyWall = "down")[ ; determine whether or not there are enough patches to support vents

    if(numSupplyVents > grid-width)[ ; If there are too many vents.

      show "numSupplyVents exceeds the number of patches on the designated wall. Simulation aborted." ; tell users why the simulation will not move forward.
      stop ; abort the simulation.

    ]
  ]

    if(ventilReturnWall = "east" OR ventilReturnWall = "west" OR ventilReturnWall = "left" OR ventilReturnWall = "right")[ ; determine whether or not there are enough patches to support vents

    if(numReturnVents > grid-height)[ ; If there are too many vents.

      show "numReturnVents exceeds the number of patches on the designated wall. Simulation aborted." ; tell users why the simulation will not move forward.
      stop ; abort the simulation.

    ]
  ]

  if(ventilSupplyWall = "east" OR ventilSupplyWall = "west" OR ventilSupplyWall = "left" OR ventilSupplyWall = "right")[ ; determine whether or not there are enough patches to support vents

    if(numSupplyVents > grid-height)[ ; If there are too many vents.

      show "numSupplyVents exceeds the number of patches on the designated wall. Simulation aborted." ; tell users why the simulation will not move forward.
      stop ; abort the simulation.

    ]
  ]

  ]

  clear-all ; Remove all variables from a previous simulation run.
  reset-ticks ; Set up the ticks. Note: all ticks represent 1-minute intervals.

  set-default-shape people "person" ; we want to use the "person" icon for our people turtles

  ; if we're simulating the choir case study, specifically

  if (choirSim = TRUE)[

    ; readjust initial settings to reflect choir activity (these are just the mandatory initial conditions). Others exist but can be more flexible
    set num-cohorts 4
    set cohort-dur 40 ; the first cohort duration will be 40 minutes.
    set grid-height 10
    set grid-width 18
    ;set n 61
    set n_infectious 1
    set numSupplyVents 3
    set numReturnVents 1
    set rearrange-cohort TRUE
    set personPerPatch-cap 2

  ]

  ; Define the initial value of global variables.

  set totalInfected 0 ; As the system is set up, there are no infections.
  set firstInfectTime "NA" ; As the system is set up, there are no infections.
  set lastInfectTime  "NA" ; As the system is set up, there are no infections.
  set cohort 1 ; Define the value of the first cohort.
  set personDist.list (list) ; Set up the average person distance list. Initially it is an empty list.
  set totalInfected.list (list) ; Set up the list of infected susceptibles. Initially it is an empty list.
  set patchContamination.list (list) ; Set up the list of contaminated patch counts. Initially it is an empty list.
  set patchInfectiousness.list (list) ; Set up the list of patch infectiousness. Initially it is an empty list.
  set cohort-endTime cohort-dur ; the first endTime value will just be cohort-dur.
  set num_symptomatic 0 ; As the system is set up, there are no symptomatic individuals
  set num_asymptomatic 0 ; As the system is set up, there are no asymptomatic individuals

  ; Set up parameters for lognormal distributions based on instructions described in Exercise 7 of Ch. 15 of Agent-based and Individual-based Modeling: A Practical Introduction (Railsback & Grimm 2011).

  set beta_symp ln (1 + ((cough_spread-dist.sd ^ 2)/(cough_spread-dist.mean ^ 2))) ;This code defines the beta to be used in the step function (to draw from a lognormal distribution)
  set M_symp (ln(cough_spread-dist.mean) - (beta_symp / 2)) ;This code defines the M to be used in the step function (to draw from a lognormal distribution)
  set S_symp sqrt beta_symp ;This code defines the S to be used in the step function (to draw from a lognormal distribution)

  set beta_asymp ln (1 + ((speak_spread-dist.sd ^ 2)/(speak_spread-dist.mean ^ 2))) ;This code defines the beta to be used in the step function (to draw from a lognormal distribution)
  set M_asymp (ln(speak_spread-dist.mean) - (beta_asymp / 2)) ;This code defines the M to be used in the step function (to draw from a lognormal distribution)
  set S_asymp sqrt beta_asymp ;This code defines the S to be used in the step function (to draw from a lognormal distribution)

  set beta_speak_droplets ln (1 + ((speak_dropletNum.sd ^ 2)/(speak_dropletNum.mean ^ 2))) ;This code defines the beta to be used in the step function (to draw from a lognormal distribution)
  set M_speak_droplets (ln(speak_dropletNum.mean) - (beta_speak_droplets / 2)) ;This code defines the M to be used in the step function (to draw from a lognormal distribution)
  set S_speak_droplets sqrt beta_speak_droplets ;This code defines the S to be used in the step function (to draw from a lognormal distribution)

  set beta_cough_droplets ln (1 + ((cough_dropletNum.sd ^ 2)/(cough_dropletNum.mean ^ 2))) ;This code defines the beta to be used in the step function (to draw from a lognormal distribution)
  set M_cough_droplets (ln(cough_dropletNum.mean) - (beta_cough_droplets / 2)) ;This code defines the M to be used in the step function (to draw from a lognormal distribution)
  set S_cough_droplets sqrt beta_cough_droplets ;This code defines the S to be used in the step function (to draw from a lognormal distribution)

  ; Define terminal downward velocities for respiratory droplets with diameters ranging from 3 micrometers to 750 micrometers.

  set Vt_diam3 0.000271406 ; The terminal velocity (m/s) of a respiratory droplet with a 3-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam6 0.001085622 ; The terminal velocity (m/s) of a respiratory droplet with a 6-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam12 0.004342489 ; The terminal velocity (m/s) of a respiratory droplet with a 12-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam20 0.012062469 ; The terminal velocity (m/s) of a respiratory droplet with a 20-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam28 0.02364244 ; The terminal velocity (m/s) of a respiratory droplet with a 28-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam36 0.0390824 ; The terminal velocity (m/s) of a respiratory droplet with a 36-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam45 0.06106625 ; The terminal velocity (m/s) of a respiratory droplet with a 45-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam62.5 0.11779755 ; The terminal velocity (m/s) of a respiratory droplet with a 62.5-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam87.5 0.230883198 ; The terminal velocity (m/s) of a respiratory droplet with a 87.5-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam112.5 0.381664063 ; The terminal velocity (m/s) of a respiratory droplet with a 112.5-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam137.5 0.570140143 ; The terminal velocity (m/s) of a respiratory droplet with a 137.5-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam175 0.923532793 ; The terminal velocity (m/s) of a respiratory droplet with a 175-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam225 1.52665625 ; The terminal velocity (m/s) of a respiratory droplet with a 225-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam375 4.240711806 ; The terminal velocity (m/s) of a respiratory droplet with a 375-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).
  set Vt_diam750 16.96284722 ; The terminal velocity (m/s) of a respiratory droplet with a 750-micrometer diameter, calculated from equations given by Anchordopqui & Chudnovsky (2020) (Preprint available at https://arxiv.org/pdf/2003.13689.pdf).


  ; Define the grid cells

  resize-world 0 (grid-width - 1) 0 (grid-height - 1) ;setup the grid (world). Each patch is a grid cell. Note: the minus ones here are because we must include zeroes on the grid.
  set-patch-size 40 ; make the simulation view bigger

  ask patches [

    set can-sprout? TRUE ; at the beginning of the simulation all patches can sprout people.
    set person-count 0 ; at the beginning of the simulation no people have sprouted.
    ;set person-dist 0 ; at the beginning of the simulation no people have sprouted. ;(IMPROVEMENT NOTE:: This variable may be added to facilitate movement while maintaining social distancing)
    set pcolor 89 ; all non-contaminated patches are this shade of cyan
    set returnVent FALSE ; all patches initially start out as a non-vent patch.
    set supplyVent FALSE ; all patches initially start out as a non-vent patch.

  ]

  ; Define the vent locations (only if airflow = TRUE)
  ;; vents will be randomly distributed on the designated wall

  if(ventilation = TRUE)[

    ifelse (equallySpaceVents = TRUE)

    [ ; if vents should be equally-spaced along walls (as much as possible, anyway).

  if(ventilReturnWall = "north" OR ventilReturnWall = "up" )[ ; determine where the vents are located

       let returnCount 0 ; create a variable that will be updated for each vent (keeps track of what vent is being placed).

       repeat numReturnVents [ ; place one vent at a time

        set returnCount returnCount + 1 ; update the vent-tracking variable

        let return.x (((grid-width / numReturnVents) * returnCount ) + ((grid-width / numReturnVents) * (returnCount - 1))) / 2 ; calculate the x coordinate for the vent.

        ask patches with [pxcor = (min-pxcor + floor return.x) AND pycor = max-pycor][ ; create the vent at the appropriate space.

          set returnVent TRUE

        ]

      ]

  ]

  if(ventilReturnWall = "south" OR ventilReturnWall = "down" )[ ; determine where the vents are located

       let returnCount 0 ; create a variable that will be updated for each vent (keeps track of what vent is being placed).

       repeat numReturnVents [ ; place one vent at a time

        set returnCount returnCount + 1 ; update the vent-tracking variable

        let return.x (((grid-width / numReturnVents) * returnCount ) + ((grid-width / numReturnVents) * (returnCount - 1))) / 2 ; calculate the x coordinate for the vent.

        ask patches with [pxcor = (min-pxcor + floor return.x) AND pycor = min-pycor][ ; create the vent at the appropriate space.

          set returnVent TRUE

        ]

      ]

  ]

  if(ventilReturnWall = "east" OR ventilReturnWall = "right" )[ ; determine where the vents are located

       let returnCount 0 ; create a variable that will be updated for each vent (keeps track of what vent is being placed).

       repeat numReturnVents [ ; place one vent at a time

        set returnCount returnCount + 1 ; update the vent-tracking variable

        let return.y (((grid-height / numReturnVents) * returnCount ) + ((grid-height / numReturnVents) * (returnCount - 1))) / 2 ; calculate the y coordinate for the vent.

        ask patches with [pycor = (min-pycor + floor return.y) AND pxcor = max-pxcor][ ; create the vent at the appropriate space.

          set returnVent TRUE

        ]

      ]

  ]

  if(ventilReturnWall = "west" OR ventilReturnWall = "left" )[ ; determine where the vents are located

       let returnCount 0 ; create a variable that will be updated for each vent (keeps track of what vent is being placed).

       repeat numReturnVents [ ; place one vent at a time

        set returnCount returnCount + 1 ; update the vent-tracking variable

        let return.y (((grid-height / numReturnVents) * returnCount ) + ((grid-height / numReturnVents) * (returnCount - 1))) / 2 ; calculate the y coordinate for the vent.

        ask patches with [pycor = (min-pycor + floor return.y) AND pxcor = min-pxcor][ ; create the vent at the appropriate space.

          set returnVent TRUE

        ]

      ]

  ]

    if(ventilSupplyWall = "north" OR ventilSupplyWall = "up" )[ ; determine where the vents are located

       let supplyCount 0 ; create a variable that will be updated for each vent (keeps track of what vent is being placed).

       repeat numSupplyVents [ ; place one vent at a time

        set supplyCount supplyCount + 1 ; update the vent-tracking variable

        let supply.x (((grid-width / numSupplyVents) * supplyCount ) + ((grid-width / numSupplyVents) * (supplyCount - 1))) / 2 ; calculate the x coordinate for the vent.

        ask patches with [pxcor = (min-pxcor + floor supply.x) AND pycor = max-pycor][ ; create the vent at the appropriate space.

          set supplyVent TRUE

        ]

      ]

  ]

     if(ventilSupplyWall = "south" OR ventilSupplyWall = "down" )[ ; determine where the vents are located

       let supplyCount 0 ; create a variable that will be updated for each vent (keeps track of what vent is being placed).

       repeat numSupplyVents [ ; place one vent at a time

        set supplyCount supplyCount + 1 ; update the vent-tracking variable

        let supply.x (((grid-width / numSupplyVents) * supplyCount ) + ((grid-width / numSupplyVents) * (supplyCount - 1))) / 2 ; calculate the x coordinate for the vent.

        ask patches with [pxcor = (min-pxcor + floor supply.x) AND pycor = min-pycor][ ; create the vent at the appropriate space.

          set supplyVent TRUE

        ]

      ]

  ]

  if(ventilSupplyWall = "east" OR ventilSupplyWall = "right" )[ ; determine where the vents are located

     let supplyCount 0 ; create a variable that will be updated for each vent (keeps track of what vent is being placed).

       repeat numSupplyVents [ ; place one vent at a time

        set supplyCount supplyCount + 1 ; update the vent-tracking variable

        let supply.y (((grid-height / numSupplyVents) * supplyCount ) + ((grid-height / numSupplyVents) * (supplyCount - 1))) / 2 ; calculate the y coordinate for the vent.

        ask patches with [pycor = (min-pycor + floor supply.y) AND pxcor = max-pxcor][ ; create the vent at the appropriate space.

          set supplyVent TRUE

        ]

      ]

  ]

  if(ventilSupplyWall = "west" OR ventilSupplyWall = "left" )[ ; determine where the vents are located

     let supplyCount 0 ; create a variable that will be updated for each vent (keeps track of what vent is being placed).

       repeat numSupplyVents [ ; place one vent at a time

        set supplyCount supplyCount + 1 ; update the vent-tracking variable

        let supply.y (((grid-height / numSupplyVents) * supplyCount ) + ((grid-height / numSupplyVents) * (supplyCount - 1))) / 2 ; calculate the y coordinate for the vent.

        ask patches with [pycor = (min-pycor + floor supply.y) AND pxcor = min-pxcor][ ; create the vent at the appropriate space.

          set supplyVent TRUE

        ]

      ]

  ]


    ]


    [ ; if vents should be randomly placed along walls.

  if(ventilReturnWall = "north" OR ventilReturnWall = "up" )[ ; determine where the vents are located

    ask n-of numReturnVents patches with [pycor = max-pycor][ ; ask the pre-determined number of patches at the top of the world to become vents

      set returnVent TRUE

    ]

  ]

  if(ventilReturnWall = "south" OR ventilReturnWall = "down" )[ ; determine where the vents are located

    ask n-of numReturnVents patches with [pycor = min-pycor][ ; ask the pre-determined number of patches at the top of the world to become vents

      set returnVent TRUE

    ]

  ]

  if(ventilReturnWall = "east" OR ventilReturnWall = "right" )[ ; determine where the vents are located

    ask n-of numReturnVents patches with [pycor = max-pxcor][ ; ask the pre-determined number of patches at the top of the world to become vents

      set returnVent TRUE

    ]

  ]

  if(ventilReturnWall = "west" OR ventilReturnWall = "left" )[ ; determine where the vents are located

    ask n-of numReturnVents patches with [pycor = min-pxcor][ ; ask the pre-determined number of patches at the top of the world to become vents

      set returnVent TRUE

    ]

  ]

    if(ventilSupplyWall = "north" OR ventilSupplyWall = "up" )[ ; determine where the vents are located

    ask n-of numSupplyVents patches with [pycor = max-pycor][ ; ask the pre-determined number of patches at the top of the world to become vents

      set supplyVent TRUE

    ]

  ]

  if(ventilSupplyWall = "south" OR ventilSupplyWall = "down" )[ ; determine where the vents are located

    ask n-of numSupplyVents patches with [pycor = min-pycor][ ; ask the pre-determined number of patches at the top of the world to become vents

      set supplyVent TRUE

    ]

  ]

  if(ventilSupplyWall = "east" OR ventilSupplyWall = "right" )[ ; determine where the vents are located

    ask n-of numSupplyVents patches with [pycor = max-pxcor][ ; ask the pre-determined number of patches at the top of the world to become vents

      set supplyVent TRUE

    ]

  ]

  if(ventilSupplyWall = "west" OR ventilSupplyWall = "left" )[ ; determine where the vents are located

    ask n-of numSupplyVents patches with [pycor = min-pxcor][ ; ask the pre-determined number of patches at the top of the world to become vents

      set supplyVent TRUE

    ]

  ]

    ]

  ask patches with [returnVent = FALSE][

   sprout-airArrows 1 ; all non-returnVent patches should sprout an airArrow

  ]


  ask airArrows [ ; we ask all the airArrows to determine the closest returnVent patch and set the heading towards it.

    ifelse (showArrows = FALSE)
      [  ;arrows will only be hidden if we want them to be.
       set hidden? TRUE ; we make it so that the arrows don't show up on the display
      ]
      [; if we can see them
        set color 97 ; most arrows will be this shade of blue
        if( [supplyVent] of patch-here = TRUE)[ ; if arrows are on a supply vent, they'll be a distinct shade.

          set color 105 ; change these arrows to this shade of blue

        ]

      ]

    if (count patches with [returnVent = TRUE] > 0)[ ; this prevents an error from occuring if numReturnVents = 0

    let distList0 (list) ; create empty list to hold vent distances

    ask patches with [returnVent = TRUE] [ ; ask returnVent patches to calculate the distance to the airArrow.

       set distList0 lput (distance myself) distList0 ; add distances to distList

     ]

    let faceMe.x 0 ; create empty object to hold x coordinate of the closest returnVent patch
    let faceMe.y 0 ; create empty object to hold y coordinate of the closest returnVent patch

    ask one-of patches with [returnVent = TRUE] in-radius min distList0 [ ; ask the closest returnVent patch to provide x and y information.  If there is > 1 returnVent patch within the defined radius, we randomly select one.

        set faceMe.x pxcor ; update x coordinate information
        set faceMe.y pycor ; update y coordinate information

      ]

    facexy faceMe.x faceMe.y ;set heading to towards the closest returnVent. If there is > 1 returnVent patch within the defined radius, we randomly select one to set the heading towards.

    ]

  ]


  ]

  ; sprout the people


  if(n > 0 AND n > n_infectious)[ ; throw this if statement here to prevent an error when n = 0 or when there are supposed to be more infectious people than total ones. In this case, the simulation will abort when it tries to go (without triggering a fatal error)

    materialize ; call the materialize action to have the people spawn at designated patches

    ; define person variables and start the infection

    ask people [

      ;I'm not quite sure how I want people to face patches the may infect. I've made t toggleable with the face-northward button.
      ifelse (face-northward = TRUE) ; determine if you want people to look ahead (within a  range of 90 degrees) or face a random direction.
        [set heading 360 + ((random-float 90) - 45)] ; Ensure the person is forward-facing (i.e., between 315 and 45 degrees) (as if they were watching someone/something at the front of the room (e.g., a professor, movie, etc.)
        [set heading random-float 360] ; have all turtles initialize facing a random direction

      set color 64 ; make all people green
      set cohort-person cohort ; denote what cohort this person belongs to.
      set infected? FALSE ; initially sprout only healthy people
      set exposure-risk 1 ; individuals default to complete susceptibility to infectious agents.
      set expectorate-risk 1 ; individuals default to an expectorate-risk of 1
      set infectious? FALSE ; initially sprout only healthy people
      set mask? FALSE ; initially sprout people not wearing masks
      set symptomatic? FALSE ; initially sprout only healthy people
      set timesteps_exposed 0 ; initially there is no chance of being on a contaminated patch.

      let risk-change? random-float 1 ; pull a random number between zero and 1

      ;update risk of exposure for susceptible agents

      if(mod_group = "sus" OR mod_group = "sus_inf")[ ;if susceptibles are designated as a group to have their risk updated

        if(risk-change? <= mod-proportion)[ ; if risk-change is less than or equal to mod-proportion, then peoples exposure-risk is changed to reflect the maskRisk-mod value

          set exposure-risk maskRisk-mod ; update exposure-risk ; Note that this is the same as multiplying exposure-risk by maskRisk-mod because exposure-risk defaults to 1.
          set mask? TRUE ; individuals' expectorate/exposure probability only changes if they're wearing a mask.

        ]
      ]

      ;update the risk of expectorating, if applicable.

      if(mod_group = "inf" OR mod_group = "sus_inf")[ ;if infectious agents are designated as a group to have their risk updated

        if(risk-change? <= mod-proportion)[ ; if risk-change is less than or equal to mod-proportion, then peoples expectorte-risk is changed to reflect the maskRisk-mod * 1 value

          set expectorate-risk maskRisk-mod ; update exposure-risk ; Note that this is the same as multiplying exposure-risk by maskRisk-mod because exposure-risk defaults to 1.
          set mask? TRUE ; individuals' expectorate/exposure probability only changes if they're wearing a mask.

        ]

      ]

    ]

    ; designate infectious agents

    ask n-of n_infectious people with [infectious? = FALSE][  ; we randomly choose n_infectious people from the group of susceptibles to be the infectious individual(s) and spread the infection

      set color 135 ;asymptomatic infectious people will be this shade of pink
      set infectious? TRUE ; obviously infectious individuals will have this variable change to TRUE. Note that for these people, we leave symptomatic? FALSE

    ]

    ask people with [infectious? = TRUE][  ; infectious people have a symp-pr chance to become symptomatic.

      let becomeSymptomatic random-float 1 ; pull a random number between zero and 1

      if(becomeSymptomatic <= symp-pr)[ ; if becomeSymptomatic is less than or equal to symp-pr, we designate them as "symptomatic"

        set color 15 ; symptomatic individuals will be this shade of red
        set symptomatic? TRUE ; these individuals will be symptomatic

      ]
    ]

    ;update global variables now that agents have spawned.

    set num_completelySusceptible count people with [infectious? = FALSE AND exposure-risk = 1] ; count the number of completely susceptible individuals (who are not infectious) in the cohort
    set num_reducedSusceptible count people with [infectious? = FALSE AND exposure-risk < 1] ; count the number of individuals who are not infectious and have reduced susceptibility in the cohort
    set num_asymptomatic count people with [infectious? = TRUE AND symptomatic? = FALSE] ; count the number of individuals that are infectious, but completely asymptomatic
    set num_symptomatic count people with [infectious? = TRUE AND symptomatic? = TRUE] ; count the number of individuals that are infectious, but may cough/sneeze during the simulation
  ]

  if(n = 0)[ ; If there are 0 individuals present, the simulation will abort

    show "No individuals present in simulation. Simulation aborted." ; tell users why the simulation will not move forward.
    stop ; abort the simulation.

  ]

  if(n <= n_infectious)[ ; there needs to be at least one sussceptible person in the model, otherwise the simulation will abort

    show "No susceptible individuals present in simulation. Simulation aborted." ; tell users why the simulation will not move forward.
    stop ; abort the simulation.

  ]

  if((num_asymptomatic + num_symptomatic) = 0)[ ; there needs to be at least one infectious person in the first cohort of the model, otherwise the simulation will abort

    show "No infectious individuals present in simulation. Simulation aborted." ; tell users why the simulation will not move forward.
    stop ; abort the simulation.

  ]

  assess-distance ; call the assess-distance action to add a new observation to personDist.list ; since people are unmoving in this model, we only need to assess distances one time.

end

to go

  if (totalInfected = (num_completelySusceptible + num_reducedSusceptible) AND cohort = num-cohorts)[ ; if all susceptible people have been infected and there won't be any more added, then there's no need to continue the simulation

    set lastInfectTime ticks ; return the time that all individuals were infected
    show "No susceptible individuals present in simulation. Simulation aborted." ; tell users why the simulation will not move forward.
    stop ; abort the simulation.

  ]

  if (evalQuanta = TRUE) [ ;if people just want to test the proportion of individuals infected if droplets were introduced in a well-mixed area, then that's what we do.

    testQuanta
    stop

  ]

  removeDroplets ; remove a proportion of the droplets that accumulated during the previous tick

  if (ventilation = TRUE) [ ;if there is airflow in the simulation, we call the move-air action

    move-air

  ]

  diffuseDroplets ; diffuse droplets to nearby patches

  color-patches ; update patch color

  tick ; make the tick counter go up. Remember, each tick represents one minute passing

  infect ; call the move action to spread infectious droplets in patches and potentially infect susceptible people

  color-patches ; update patch color

  ;assess-distance ; call the assess-distance action to add a new observation to personDist.list.

  ;set personDist.list lput avg.dist personDist.list; keep track of the mean inter-personal distance value at each time step.

  assess-contamination ; call on the assess-contamination action to update the avg.PatchInfectiousness value at each time step.

  assess-exposures ; call the assess-exposures action to update the avg.exposureTime variable.

  set totalInfected.list lput totalInfected totalInfected.list; keep track of the totalInfected value at each time step.

  set patchInfectiousness.list lput avg.PatchInfectiousness patchInfectiousness.list; keep track of the mean patch infectiousness value at each time step. Note that this is assessed prior to removing droplets in patches because we're interested in the risk that individuals experienced.

  set patchContamination.list lput (count patches with [totalDroplets > 0]) patchContamination.list; keep track of the number of contaminated patches at each time step.Note that this is assessed prior to removing droplets in patches because we're interested in the risk that individuals experienced.

  if(choirSim = TRUE)

  [ ; if we're explicitly simulating the choir case study

    if(ticks = 40 AND num-cohorts > 1)[ ; if we reach the number of minutes for the time spent in the first practice room, we adjust cohort-dur to reflect time spent in the second rooms.

      ask people [ ; ask people to remember the first patch the started on so that they could return to it later
        set rememberPatch patch-here
      ]
      set cohort-dur 50 ; reset cohort-dur

    ]

    if(ticks = (40 + 50) AND num-cohorts > 2)[ ; if we reach the maximum number of minutes for the time spent in the first AND second practice rooms, we adjust cohort-dur to reflect time spent during break.

      set cohort-dur 15 ; reset cohort-dur

    ]

    if(ticks = (40 + 50 + 15) AND num-cohorts > 3)[ ; if we reach the maximum number of minutes for the time spent in the first AND second practice rooms AND the break, we adjust cohort-dur to reflect time spent in the final location (which is the same location as the first one).

      set cohort-dur 45 ; reset cohort-dur

    ]

  ]

   if(cohort >= num-cohorts AND ticks >= cohort-endTime)[ ; if we reach the maximum number of cohorts AND the maximum cohort-endTime, we stop the simulation.

      stop
    ]

    if(ticks >= cohort-endTime)[ ;if we reach the set cohort-endTime, but NOT the maximum number of cohorts, we replace the people

      cohort-replace ; call the cohort-replace action

    ]

    if(choirSim = TRUE AND ticks = (40 + 50 + 15) AND num-cohorts >= 4 )[ ; if we're explicitly simulating the choir case study and it's the final movement time point

    ask people [move-to rememberPatch] ; have people move back to their original patch

  ]


end

to assess-distance ; we want to keep track of the per-capita average inter-personal distances over the course of the simulation

  let distList (list) ; create a temporary list to hold distances to each person on a given time step

  ask people [ ; we ask all the people
    ask other people [ ; to ask all other people to calculate the distance to the originally-asked person.

      set distList lput (distance myself) distList ; add distances to distList

    ]
  ]

  set avg.dist mean distList ; calculate the mean
  set avg.dist precision avg.dist 2 ; round to decimal points
  ;set personDist.list lput avg.dist personDist.list ; add the mean value of the averaged interpersonal distance values to personDist.list

end

to assess-contamination ; we want to keep track of the average infectiousness of patchs in the simulation.

  let infectiousnessList (list) ; create a temporary list to hold total infectiousness count values on a given timestep

  ask patches[

    set infectiousnessList lput (transmissionRisk) infectiousnessList ; ask each patch to add the probability of infection given exposure to virions present to this list.

  ]

  set avg.PatchInfectiousness mean infectiousnessList ; calculate the mean
  set avg.PatchInfectiousness precision avg.PatchInfectiousness 2 ; round to decimal points

end

to assess-exposures ; update the avg.exposureTime variables

  let exposureTimeCountList (list) ; create a temporary list to hold total exposure-timestep count values on a given timestep

  ask people with [infectious? = FALSE] [ ; we ask only susceptible people. Otherwise this average will be inflated as infectious people are ALWAYS within any cone of exposure that they produce.

    set exposureTimeCountList lput (timesteps_exposed) exposureTimeCountList ; ask each person to add the number of timesteps where they were exposed to this list.

  ]

  set avg.exposureTime mean exposureTimeCountList ; calculate the mean
  set avg.exposureTime precision avg.exposureTime 2 ; round to decimal points

end

to cohort-rearrange ; Instead of replacing all cohort members, we move them around the world. This is effectively a modified "materialize" action.


  repeat count people [ ; we use a repeat here to move people one person at a time to ensure that moving people observe social distancing when possible

      ifelse (count patches with [can-sprout? = TRUE] > 0) ; we first ask if there are any patches that can hold people while maintaining social distancing
       [ ;if at least one viable patch DOES exist, we ask a person to move there

        ask one-of people with [cohort-person = (cohort - 1)] [ ; ask a person that hasn't moved yet.

          set cohort-person cohort ; update the cohort number (this is important because it is used below to assess distances between rearranged people and patches)
          move-to one-of patches with [can-sprout? = TRUE] ; move to an available patch

          ask patch-here [
            set person-count (person-count + 1) ; update the person-count patch variable
          ]

        ]


        ask patches [ ; we ask all the patches to calculate the distance to newly-moved people. If the minimum distance is less than social-distance, we tell them to set can-sprout false.

          let distList1 (list) ; create empty list to hold person distances

          ask people with [cohort-person = cohort] [ ; to ask all people that have already moved to calculate the distance to the originally-asked patch.

            set distList1 lput (distance myself) distList1 ; add distances to distList

          ]

          if( (min distList1) < social-distance)[ ; if the minimum distance to a person is less than social-distance, patches must turn off their ability to sprout more people (Note: if social-distance = 0, can-sprout? remains as TRUE because no distance will be < 0).

            set can-sprout? FALSE ; turn off patches' ability to sprout more people

          ]

          if(person-count = personPerPatch-cap)[ ; if the patch has reached the maximum number of people it can hold (Note: this will only be relevant here if social-distance = 0).

            set can-sprout? FALSE ; turn off patches' ability to sprout more people

          ]

        ]
       ]
     [  ;if No viable patches exist, we ask people to spawn on the one of the least-populated patches

        let popList (list) ; create an empty list that we will append turtle counts to

        ask patches [

          set popList lput person-count popList ; ask all patches to append the count of people on them to the list

        ]

       ask one-of people with [cohort-person = (cohort - 1)] [ ; ask a person that hasn't moved yet.

          set cohort-person cohort ; update the cohort number (this is important because it is used below to assess distances between rearranged people and patches)
          move-to one-of patches with [person-count = min popList] ; move to an available patch with the lowest person-count values (because we assume turtles don't want to cluster)

          ask patch-here [
            set person-count (person-count + 1) ; update the person-count patch variable
          ]

    ]
  ]
  ]

end

to cohort-replace ; The cohort-replace action is essentially the setup action, without redefining certain global and local variables. Note that all indivivduals in cohorts after the first one will be susceptible (unless rearrange-cohort = TRUE, in which case cohort-replace effectively becomes a rough proxy for movement of people within the room).

  set cohort (cohort + 1) ; Add 1 to the cohort value.
  set cohort-endTime (cohort-endTime + cohort-dur) ; Add to the cohort-endTime value.

  ask patches [

    set can-sprout? TRUE ; at the beginning of the simulation all patches can sprout people.
    set person-count 0 ; at the beginning of the simulation no people have sprouted.
    ;set person-dist 0 ; at the beginning of the simulation no people have sprouted. ;(IMPROVEMENT NOTE:: This variable may be added to facilitate movement while maintaining social distancing)

  ]

  ; sprout the people

  ifelse (rearrange-cohort = TRUE)

   [ cohort-rearrange ] ; if TRUE, we just call the cohort-rearrange action

   [ ; if FALSE, we spawn completely new turtles

    clear-turtles ; remove all turtles from the previous cohort
    materialize ; call the materialize action to have the people spawn at designated patches

    ; define person variables and start the infection

    ask people [

      set color 64 ; make all people green
      set cohort-person cohort ; denote what cohort this person belongs to.
      set infected? FALSE ; initially sprout only healthy people
      set exposure-risk 1 ; individuals default to complete susceptibility to infectious agents.
      set expectorate-risk 1 ; individuals default to an expectorate-risk of 1
      set infectious? FALSE ; initially sprout only healthy people
      set mask? FALSE ; initially sprout people not wearing masks
      set symptomatic? FALSE ; initially sprout only healthy people
      set timesteps_exposed 0 ; initially there is no chance of being on a contaminated patch.

      let risk-change? random-float 1 ; pull a random number between zero and 1

      ;update risk of exposure for susceptible agents

      if(mod_group = "sus" OR mod_group = "sus_inf")[ ;if susceptibles are designated as a group to have their risk updated

        if(risk-change? <= mod-proportion)[ ; if risk-change is less than or equal to mod-proportion, then peoples exposure-risk is changed to reflect the maskRisk-mod value

          set exposure-risk maskRisk-mod ; update exposure-risk ; Note that this is the same as multiplying exposure-risk by maskRisk-mod because exposure-risk defaults to 1.
          set mask? TRUE ; individuals' expectorate/exposure probability only changes if they're wearing a mask.

        ]
      ]

    ]

    ;update global variables now that agents have spawned.

    set num_completelySusceptible count people with [infectious? = FALSE AND exposure-risk = 1] ; count the number of completely susceptible individuals (who are not infectious) in the cohort
    set num_reducedSusceptible count people with [infectious? = FALSE AND exposure-risk < 1] ; count the number of individuals who are not infectious and have reduced susceptibility in the cohort
    set num_asymptomatic count people with [infectious? = TRUE AND symptomatic? = FALSE] ; count the number of individuals that are infectious, but completely asymptomatic
    set num_symptomatic count people with [infectious? = TRUE AND symptomatic? = TRUE] ; count the number of individuals that are infectious, but may cough/sneeze during the simulation

]

    ask people [ ; update headings

      ;I'm not quite sure how I want people to face patches the may infect. I've made t toggleable with the face-northward button.
      ifelse (face-northward = TRUE) ; determine if you want people to look ahead (within a  range of 90 degrees) or face a random direction.
        [set heading 360 + ((random-float 90) - 45)] ; Ensure the person is forward-facing (i.e., between 315 and 45 degrees) (as if they were watching someone/something at the front of the room (e.g., a professor, movie, etc.)
        [set heading random-float 360] ; have all turtles initialize facing a random direction

  ]

end

to color-patches ; for viewing purposes, we want to scale patch coloration based on infection risk

  ask patches [
    ifelse (transmissionRisk > 0)
    [set pcolor scale-color magenta (transmissionRisk * (vol_B / (1 * 1 * expectorate-height))) 0.99 0] ; all contaminated patches will be a shade of pink. The darker the color, the more infectious the patch is.
    [set pcolor 89] ; all non-contaminated patches are this shade of cyan
  ]

end

to diffuseDroplets ; diffuse droplet counts across to surrounding patches according to the diffusionRate (i.e., rate at which droplets of all sizes diffuse out of patches).

  ask patches [ ; ensure additionalDroplet variables are cleared

    set additionalDroplets_size3  0  ; reset the additional droplet counter
    set additionalDroplets_size6   0  ; reset the additional droplet counter
    set additionalDroplets_size12  0  ; reset the additional droplet counter
    set additionalDroplets_size20  0  ; reset the additional droplet counter
    set additionalDroplets_size28  0  ; reset the additional droplet counter
    set additionalDroplets_size36  0  ; reset the additional droplet counter
    set additionalDroplets_size45  0  ; reset the additional droplet counter
    set additionalDroplets_size62.5  0  ; reset the additional droplet counter
    set additionalDroplets_size87.5  0  ; reset the additional droplet counter
    set additionalDroplets_size112.5  0  ; reset the additional droplet counter
    set additionalDroplets_size137.5  0  ; reset the additional droplet counter
    set additionalDroplets_size175  0  ; reset the additional droplet counter
    set additionalDroplets_size225  0  ; reset the additional droplet counter
    set additionalDroplets_size375  0  ; reset the additional droplet counter
    set additionalDroplets_size750 0  ; reset the additional droplet counter

  ]

  ask patches [ ; counts to the additionalDroplet variables of neighbors

    let numNeighbors count neighbors ; count the number of neighbors each patch has

    ;;Note that below we use the same "ReturnCount" nomenclature as that used in the move-air action. This is purely due to (admittedly-lazy) copy-pasting from that action. Return vents are not actually directly involved here.

    let dropletSize3ReturnCount (droplets_size3 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize6ReturnCount (droplets_size6 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize12ReturnCount (droplets_size12 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize20ReturnCount (droplets_size20 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize28ReturnCount (droplets_size28 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize36ReturnCount (droplets_size36 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize45ReturnCount (droplets_size45 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize62.5ReturnCount (droplets_size62.5 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize87.5ReturnCount (droplets_size87.5 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize112.5ReturnCount (droplets_size112.5 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize137.5ReturnCount (droplets_size137.5 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize175ReturnCount (droplets_size175 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize225ReturnCount (droplets_size225 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize375ReturnCount (droplets_size375 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.
  let dropletSize750ReturnCount (droplets_size750 * diffusionRate)  ; determine the number of droplets of this size that will be removed from this patch due to diffusion.

    ask neighbors [ ; add droplets to neighbors

    set additionalDroplets_size3  additionalDroplets_size3 + (dropletSize3ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size6  additionalDroplets_size6 + (dropletSize6ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size12 additionalDroplets_size12 + (dropletSize12ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size20 additionalDroplets_size20 + (dropletSize20ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size28 additionalDroplets_size28 + (dropletSize28ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size36 additionalDroplets_size36 + (dropletSize36ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size45 additionalDroplets_size45 + (dropletSize45ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size62.5 additionalDroplets_size62.5 + (dropletSize62.5ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size87.5 additionalDroplets_size87.5 + (dropletSize87.5ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size112.5 additionalDroplets_size112.5 + (dropletSize112.5ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size137.5 additionalDroplets_size137.5 + (dropletSize137.5ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size175 additionalDroplets_size175 + (dropletSize175ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size225 additionalDroplets_size225 + (dropletSize225ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size375 additionalDroplets_size375 + (dropletSize375ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors
    set additionalDroplets_size750 additionalDroplets_size750 + (dropletSize750ReturnCount / numNeighbors)   ; The number of droplets is evenly distributed among neighbors


    ]

    ; remove droplets from the current patch

    set droplets_size3 droplets_size3 - dropletSize3ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size6 droplets_size6 - dropletSize6ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size12 droplets_size12 - dropletSize12ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size20 droplets_size20 - dropletSize20ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size28 droplets_size28 - dropletSize28ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size36 droplets_size36 - dropletSize36ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size45 droplets_size45 - dropletSize45ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size62.5 droplets_size62.5 - dropletSize62.5ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size87.5 droplets_size87.5 - dropletSize87.5ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size112.5 droplets_size112.5 - dropletSize112.5ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size137.5 droplets_size137.5 - dropletSize137.5ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size175 droplets_size175 - dropletSize175ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size225 droplets_size225 - dropletSize225ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size375 droplets_size375 - dropletSize375ReturnCount ; remove the appropriate number of droplets of this size
    set droplets_size750 droplets_size750 - dropletSize750ReturnCount ; remove the appropriate number of droplets of this size

  ]


  ;Now we officially add droplets, and update total droplet counts

    ask patches [ ; update droplet counts

    set droplets_size3 droplets_size3 + additionalDroplets_size3 ; add the appropriate number of droplets of this size
    set droplets_size6 droplets_size6 + additionalDroplets_size6 ; add the appropriate number of droplets of this size
    set droplets_size12 droplets_size12 + additionalDroplets_size12 ; add the appropriate number of droplets of this size
    set droplets_size20 droplets_size20 + additionalDroplets_size20 ; add the appropriate number of droplets of this size
    set droplets_size28 droplets_size28 + additionalDroplets_size28 ; add the appropriate number of droplets of this size
    set droplets_size36 droplets_size36 + additionalDroplets_size36 ; add the appropriate number of droplets of this size
    set droplets_size45 droplets_size45 + additionalDroplets_size45 ; add the appropriate number of droplets of this size
    set droplets_size62.5 droplets_size62.5 + additionalDroplets_size62.5 ; add the appropriate number of droplets of this size
    set droplets_size87.5 droplets_size87.5 + additionalDroplets_size87.5 ; add the appropriate number of droplets of this size
    set droplets_size112.5 droplets_size112.5 + additionalDroplets_size112.5 ; add the appropriate number of droplets of this size
    set droplets_size137.5 droplets_size137.5 + additionalDroplets_size137.5 ; add the appropriate number of droplets of this size
    set droplets_size175 droplets_size175 + additionalDroplets_size175 ; add the appropriate number of droplets of this size
    set droplets_size225 droplets_size225 + additionalDroplets_size225 ; add the appropriate number of droplets of this size
    set droplets_size375 droplets_size375 + additionalDroplets_size375 ; add the appropriate number of droplets of this size
    set droplets_size750 droplets_size750 + additionalDroplets_size750 ; add the appropriate number of droplets of this size

    set additionalDroplets_size3  0  ; reset the additional droplet counter
    set additionalDroplets_size6   0  ; reset the additional droplet counter
    set additionalDroplets_size12  0  ; reset the additional droplet counter
    set additionalDroplets_size20  0  ; reset the additional droplet counter
    set additionalDroplets_size28  0  ; reset the additional droplet counter
    set additionalDroplets_size36  0  ; reset the additional droplet counter
    set additionalDroplets_size45  0  ; reset the additional droplet counter
    set additionalDroplets_size62.5  0  ; reset the additional droplet counter
    set additionalDroplets_size87.5  0  ; reset the additional droplet counter
    set additionalDroplets_size112.5  0  ; reset the additional droplet counter
    set additionalDroplets_size137.5  0  ; reset the additional droplet counter
    set additionalDroplets_size175  0  ; reset the additional droplet counter
    set additionalDroplets_size225  0  ; reset the additional droplet counter
    set additionalDroplets_size375  0  ; reset the additional droplet counter
    set additionalDroplets_size750 0  ; reset the additional droplet counter

    set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count
    set virionCount (((1.40947E-14 * 1000 * virionsPerML) * droplets_size3) + ((1.12758E-13 * 1000 * virionsPerML) * droplets_size6) + ((9.02064E-13 * 1000 * virionsPerML) * droplets_size12) + ((4.17622E-12 * 1000 * virionsPerML) * droplets_size20) + ((1.14595E-11 * 1000 * virionsPerML) * droplets_size28) + ((2.43557E-11 * 1000 * virionsPerML) * droplets_size36) + ((4.75698E-11 * 1000 * virionsPerML) * droplets_size45) + ((1.27448E-10 * 1000 * virionsPerML) * droplets_size62.5) + ((3.49718E-10 * 1000 * virionsPerML) * droplets_size87.5) + ((7.43277E-10 * 1000 * virionsPerML) * droplets_size112.5) + ((1.35707E-09 * 1000 * virionsPerML) * droplets_size137.5) + ((2.79774E-09 * 1000 * virionsPerML) * droplets_size175) + ((5.94622E-09 * 1000 * virionsPerML) * droplets_size225) + ((2.75288E-08 * 1000 * virionsPerML) * droplets_size375) + ((2.2023E-07 * 1000 * virionsPerML) * droplets_size750)) ; Multiply estimated virion counts per droplet by droplet counts (Note that the first term is in each equation is the mass or volume of each droplet size in kg or L (it's interchangeable because the weight of 1 liter of pure water at temperature 4 °C = 1) . We multiply this by 1000 to get the volume of each droplet size in mL.
    set transmissionRisk (virionCount * virionRisk) ;update risk of infection in patches.

  ]


end

to infect ;try to spread the infection

  ;First we contaminate patches (if applicable)

  ask people with[infectious? = TRUE] [

      let spread.dist 0 ; we have to create the spread.dist object before the following ifelse statement, but the value will be immediately changed therein
      let airflow-angle 0 ; we have to create the airflow-angle object before the following ifelse statement, but the value will be immediately changed therein
      let coughSizeDistr 0 ; we first state that we should NOT base the droplet-cluster size on the coughing size probability distribution (i.e., we base it on the speaking probability). If the infectious individual coughs, we update this number.
      let dropletNum 0 ; we have to create the dropletNum object before the following ifelse statement, but the value will be immediately changed therein

      ifelse(symptomatic? = TRUE) ; Here we ensure we sample from the appropriate log-normal distributions for symptomatic and asymptomatic people.
      [ ; if people are symptomatic

        let will_cough random-float 1 ; pull a number random number between 0 and 1 to determine if people will cough, or if the expectoration event will be parameterized to reflect speaking.

        ifelse(will_cough <= cough_frequency)[ ; if will_cough is less than or equal to cough_frequency, the resulting cone of exposure will be produced based on cough parameters. If not, the cone of exposure will be based on speaking parameters.

          ;if agents will cough
          set coughSizeDistr 1 ; we reset the cluster size distribution to reflect that infectious agents cough
          set airflow-angle cough_airflow-angle ; ensure we use the appropriate angle
          set spread.dist exp (random-normal M_symp S_symp) ; pull the distance (in meters) that infectious droplets will spread from a lognormal distribution following the procedure described by (Railsback & Grimm 2011)
          set dropletNum exp (random-normal M_cough_droplets S_cough_droplets) ; pull the number of droplets that will be added to patch counts from a lognormal distribution following the procedure described by (Railsback & Grimm 2011)
        ]
        [ ;if agents will not cough

          set airflow-angle speak_airflow-angle ; ensure we use the appropriate angle
          set spread.dist exp (random-normal M_asymp S_asymp) ; pull the distance (in meters) that infectious droplets will spread from a lognormal distribution following the procedure described by (Railsback & Grimm 2011)
          set dropletNum exp (random-normal M_speak_droplets S_speak_droplets) ; pull the number of droplets that will be added to patch counts from a lognormal distribution following the procedure described by (Railsback & Grimm 2011)
        ]

      ]
      [ ; if people are asymptomatic

        set airflow-angle speak_airflow-angle ; ensure we use the appropriate angle
        set spread.dist exp (random-normal M_asymp S_asymp) ; pull the distance (in meters) that infectious droplets will spread from a lognormal distribution following the procedure described by (Railsback & Grimm 2011)
        set dropletNum exp (random-normal M_speak_droplets S_speak_droplets) ; pull the number of droplets that will be added to patch counts from a lognormal distribution following the procedure described by (Railsback & Grimm 2011)
      ]

      if(mask? = TRUE)[ ; if infectious individuals are wearing a mask, we assume that any droplets produced will never extend past the patch containing the infectious individual.

        set spread.dist 0 ; the droplet-spread distance becomes 0.
        ;set dropletNum (dropletNum * 0.1) ; reduce the number of droplets by a fixed percentage.

      ]

      set dropletNum (dropletNum * expectorate-risk) ; if individuals are wearing masks, the number of droplets they expel is modified by expectorate-risk (Note that if they are not wearing masks expectorate-risk is 1)

      let patchCount count patches in-cone spread.dist airflow-angle ; count the number of patches in the cone of infection
      ;let contamPatchCount random patchCount ; Randomly decide the number of patches in the cone of infection over which the total number of droplets will be evenly distributed.
      ;set contamPatchCount contamPatchCount + 1 ; The random action calls an integer from 0 to (patchCount - 1). We want a number between 1 and patchCount.
      ;let dropletDistrNum (dropletNum / contamPatchCount) ; calculate the number of droplets each patch will receive. (Note that droplets are evenly distributed among contaminated patches)
      let dropletDistrNum (dropletNum / patchCount) ; calculate the number of droplets each patch will receive. (Note that droplets are evenly distributed among all patches within cones of infection)

      ;update droplet count numbers in contaminated patches

      ;ask n-of contamPatchCount patches in-cone spread.dist airflow-angle [ ;ask contamPatchCount randomly-selected patches within the cone of infection to update their droplet counts.
      ask patches in-cone spread.dist airflow-angle [ ;ask patches within the cone of infection to update their droplet counts.

        if(coughSizeDistr = 0)[; if infectious individual does NOT cough we sample from the speaking droplet-size probability distribution (Chao et al. 2009 - referenced below)

          ; Note that probability distribution is based on the mean droplet size range presented by Chao et al. (2009), which represent observed droplet-size measurements 60mm away from individuals' mouths.
          ; Chao, C.Y.H., Wan, M.P., Morawska, L., Johnson, G.R., Ritovski, Z.D., …, & Katoshevski, D. (2009). Characterization of expiration air jets and droplet size distributions immediately at the mouth opening. Aerosol Science 40(2009):122 – 133. https://doi.org/10.1016/j.jaerosci.2008.10.003.

          if(speak_dropletSizeDistr = "chao")[ ; if the droplet size distribution will be defined by the Chao et al. (2009) paper

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0.10526316))  ;Add the number of droplets multiplied by the probability that droplets will be 2-4 (i.e., mean of 3) micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0.36842105))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 0.15789474))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0.098398169))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0.059496568))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0.043478261))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0.022883295))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0.032036613))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0.027459954))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0.027459954))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0.009153318))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0.022883295))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0.009153318))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0.013729977))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0.00228833))  ;Add the number of droplets multiplied by the probability that droplets will be 500-1000 (i.e., mean of 750) micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(speak_dropletSizeDistr = "meanlog.1")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 1, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0.90221))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0.09769))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 1.00E-04))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(speak_dropletSizeDistr = "meanlog.2")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 2, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0.02034))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0.58528))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 0.38923))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0.0051))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 5.00E-05))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(speak_dropletSizeDistr = "meanlog.3")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 3, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0.00099))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 0.22496))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0.49652))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0.21661))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0.05001))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0.00979))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0.00111))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 1.00E-05))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(speak_dropletSizeDistr = "meanlog.4")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 4, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 3.00E-05))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0.00327))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0.03485))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0.11234))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0.23613))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0.46714))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0.12405))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0.01939))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0.00232))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0.00048))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(speak_dropletSizeDistr = "meanlog.5")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 5, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0.00012))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0.01131))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0.0818))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0.19145))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0.22957))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0.32449))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0.11992))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0.04132))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 2.00E-05))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]


        ]

       if(coughSizeDistr = 1)[; if infectious individual DOES cough we sample from the coughing droplet-size probability distribution (Chao et al. 2009 - referenced below)

          ; Note that probability distribution is based on the mean droplet size range presented by Chao et al. (2009), which represent observed droplet-size measurements 60mm away from individuals' mouths.
          ; Chao, C.Y.H., Wan, M.P., Morawska, L., Johnson, G.R., Ritovski, Z.D., …, & Katoshevski, D. (2009). Characterization of expiration air jets and droplet size distributions immediately at the mouth opening. Aerosol Science 40(2009):122 – 133. https://doi.org/10.1016/j.jaerosci.2008.10.003.

          if(cough_dropletSizeDistr = "chao")[ ; if the droplet size distribution will be defined by the Chao et al. (2009) paper

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0.091863517))  ;Add the number of droplets multiplied by the probability that droplets will be 2-4 (i.e., mean of 3) micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0.461942257))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 0.170603675))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0.073490814))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0.036745407))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0.015748031))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0.005249344))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0.023622047))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0.01312336))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0.026246719))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0.018372703))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0.015748031))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0.01312336))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0.023622047))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0.010498688))  ;Add the number of droplets multiplied by the probability that droplets will be 500-1000 (i.e., mean of 750) micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(cough_dropletSizeDistr = "meanlog.1")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 1, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0.90221))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0.09769))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 1.00E-04))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(cough_dropletSizeDistr = "meanlog.2")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 2, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0.02034))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0.58528))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 0.38923))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0.0051))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 5.00E-05))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(cough_dropletSizeDistr = "meanlog.3")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 3, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0.00099))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 0.22496))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0.49652))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0.21661))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0.05001))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0.00979))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0.00111))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 1.00E-05))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(cough_dropletSizeDistr = "meanlog.4")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 4, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 3.00E-05))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0.00327))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0.03485))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0.11234))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0.23613))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0.46714))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0.12405))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0.01939))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0.00232))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0.00048))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]

         if(cough_dropletSizeDistr = "meanlog.5")[ ; if the droplet size distribution will be defined by sampling from a lognormal distribution with meanlog droplet size of 5, sdlog of 0.3, and length of 1000000.

          set droplets_size3 (droplets_size3 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter <= 3 micrometers.
          set droplets_size6 (droplets_size6 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.
          set droplets_size36 (droplets_size36 + (dropletDistrNum * 0))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
          set droplets_size45 (droplets_size45 + (dropletDistrNum * 0.00012))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
          set droplets_size62.5 (droplets_size62.5 + (dropletDistrNum * 0.01131))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
          set droplets_size87.5 (droplets_size87.5 + (dropletDistrNum * 0.0818))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
          set droplets_size112.5 (droplets_size112.5 + (dropletDistrNum * 0.19145))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
          set droplets_size137.5 (droplets_size137.5 + (dropletDistrNum * 0.22957))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
          set droplets_size175 (droplets_size175 + (dropletDistrNum * 0.32449))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
          set droplets_size225 (droplets_size225 + (dropletDistrNum * 0.11992))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
          set droplets_size375 (droplets_size375 + (dropletDistrNum * 0.04132))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
          set droplets_size750 (droplets_size750 + (dropletDistrNum * 2.00E-05))  ;Add the number of droplets multiplied by the probability that droplets will have a mean diameter >= 750 micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count

          ]


        ]
      ]

  ]

      ;Now that we have contaminated patches, we can infect others. Infection probability is determined by the cumulative sum of virions in droplets in the patch multiplied by virionRisk (i.e., the risk of infection given exposure to a single virion).
      ; For SARS-CoV2, we assume there are 7 X 10^6 functional virions per milliliter of droplet fluid (Wolfel et al. 2020). We calculate the mass (in kg) of each droplet using the equation described by Anchordoqui & Chudnovski (2020), and convert kg to L, then to mL for virion estimation.
      ; Additionally, we assume that droplets are evenly distributed within patches, and that the air within which droplets are distubited is contained within 1 m X 1 m X expectorate-height m cubic area (i.e., patches represent 1 m X 1 m X expectorate-height m cubic space, and individuals within these patches must breathe in air from this space).

   ask patches with [totalDroplets > 0] [ ; call on infectious patches.

      set virionCount (((1.40947E-14 * 1000 * virionsPerML) * droplets_size3) + ((1.12758E-13 * 1000 * virionsPerML) * droplets_size6) + ((9.02064E-13 * 1000 * virionsPerML) * droplets_size12) + ((4.17622E-12 * 1000 * virionsPerML) * droplets_size20) + ((1.14595E-11 * 1000 * virionsPerML) * droplets_size28) + ((2.43557E-11 * 1000 * virionsPerML) * droplets_size36) + ((4.75698E-11 * 1000 * virionsPerML) * droplets_size45) + ((1.27448E-10 * 1000 * virionsPerML) * droplets_size62.5) + ((3.49718E-10 * 1000 * virionsPerML) * droplets_size87.5) + ((7.43277E-10 * 1000 * virionsPerML) * droplets_size112.5) + ((1.35707E-09 * 1000 * virionsPerML) * droplets_size137.5) + ((2.79774E-09 * 1000 * virionsPerML) * droplets_size175) + ((5.94622E-09 * 1000 * virionsPerML) * droplets_size225) + ((2.75288E-08 * 1000 * virionsPerML) * droplets_size375) + ((2.2023E-07 * 1000 * virionsPerML) * droplets_size750)) ; Multiply estimated virion counts per droplet by droplet counts (Note that the first term is in each equation is the mass or volume of each droplet size in kg or L (it's interchangeable because the weight of 1 liter of pure water at temperature 4 °C = 1) . We multiply this by 1000 to get the volume of each droplet size in mL.

      set transmissionRisk (virionCount * virionRisk) ;update risk of infection in patches.

      ask people-here [ ;ask people in infectious patches

        set timesteps_exposed (timesteps_exposed + 1) ; add 1 to the timesteps_exposed variable (i.e., counts the number of times individuals were on infectious patches)

        if(infected? = FALSE AND infectious? = FALSE)[ ; only uninfected people can become infected

           let will_infect? random-float 1 ; pull a random number between zero and 1
           if(will_infect? <= (exposure-risk * ([transmissionRisk] of patch-here)  * (vol_B / (1 * 1 * expectorate-height)) ))[ ; if will_infect? does not exceed the risk of infection (i.e., the probability that a person is exposed to virions in the patch AND that they breathe it in AND that inhilation will lead to infection), the person is infected.

                set color 25 ; infected people will be this shade of orange
                set infected? TRUE ; update infection status
                set totalInfected (totalInfected + 1) ; add to the totalInfected count
                set  droplets_size3AtInf  [droplets_size3] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size6AtInf  [droplets_size6] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size12AtInf  [droplets_size12] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size20AtInf  [droplets_size20] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size28AtInf  [droplets_size28] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size36AtInf  [droplets_size36] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size45AtInf  [droplets_size45] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size62.5AtInf [droplets_size62.5] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size87.5AtInf [droplets_size87.5] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size112.5AtInf [droplets_size112.5] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size137.5AtInf  [droplets_size137.5] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size175AtInf  [droplets_size175] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size225AtInf  [droplets_size225] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size375AtInf  [droplets_size375] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.
                set  droplets_size750AtInf  [droplets_size750] of patch-here ; Record the number of droplets of this size on the patch where this individual was infected.


              if(totalInfected = 1)[ ; when the totalInfected count first becomes 1, we record the time at which it happens

               set firstInfectTime ticks ; record the time of the first infection event.

              ]


            ]
            ]
    ]
    ]



end

to materialize ; spawn people directly on randomly-decided patches

  repeat n [ ; we use a repeat here to sprout people one person at a time to ensure that sprouting people observe social distancing when possible

      ifelse (count patches with [can-sprout? = TRUE] > 0) ; we first ask if there are any patches that can sprout people while maintaining social distancing
       [ ;if at least one viable patch DOES exist, we ask one of them to sprout a person

              ask one-of patches with [can-sprout? = TRUE] [ ; identify and ask a single patch capable of sprouting outside a specified distance to sprout a person

                  sprout-people 1 ; sprout a single person
                  set person-count (person-count + 1) ; update the person-count patch variable

              ]

        ask patches [ ; we ask all the patches to calculate the distance to all turtles. If the minimum distance is less than social-distance, we tell them to set can-sprout false.

          let distList1 (list) ; create empty list to hold person distances

          ask people [ ; to ask all people to calculate the distance to the originally-asked patch.

            set distList1 lput (distance myself) distList1 ; add distances to distList

          ]

          if( (min distList1) < social-distance)[ ; if the minimum distance to a person is less than social-distance, patches must turn off their ability to sprout more people (Note: if social-distance = 0, can-sprout? remains as TRUE because no distance will be < 0).

            set can-sprout? FALSE ; turn off patches' ability to sprout more people

          ]

          if(person-count = personPerPatch-cap)[ ; if the patch has reached the maximum number of people it can hold (Note: this will only be relevant here if social-distance = 0).

            set can-sprout? FALSE ; turn off patches' ability to sprout more people

          ]

        ]
       ]
     [  ;if No viable patches exist, we ask people to spawn on the one of the least-populated patches

        let popList (list) ; create an empty list that we will append turtle counts to

        ask patches [

          set popList lput person-count popList ; ask all patches to append the count of people on them to the list

        ]

        ask one-of patches with [person-count = min popList ][ ; ask one random patch with the lowest person-count values (because we assume turtles don't want to cluster) to sprout a person

           sprout-people 1 ; sprout a single person
           set person-count (person-count + 1) ; update the person-count patch variable

        ]

    ]
  ]

end

to move-air ; this is the action that moves air (and droplets suspended therein) throughout the simulation

  ;; We assume a very simple movement within an enclosed room where air moves only towards the return vent(s) and return vents return number droplets * ventil_movementRate * (1 - ventil_removalRate) droplets to supply vents.
  ;; Droplets returned to supply vents are transferred there immediately, which may not necessarily be reflective of real-world conditions.

  if (count patches with [returnVent = TRUE] > 0)[ ; this prevents an error from occuring if numReturnVents = 0

  let dropletSize3ReturnCount (sum ([droplets_size3] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize6ReturnCount (sum ([droplets_size6] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize12ReturnCount (sum ([droplets_size12] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize20ReturnCount (sum ([droplets_size20] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize28ReturnCount (sum ([droplets_size28] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize36ReturnCount (sum ([droplets_size36] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize45ReturnCount (sum ([droplets_size45] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize62.5ReturnCount (sum ([droplets_size62.5] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize87.5ReturnCount (sum ([droplets_size87.5] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize112.5ReturnCount (sum ([droplets_size112.5] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize137.5ReturnCount (sum ([droplets_size137.5] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize175ReturnCount (sum ([droplets_size175] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize225ReturnCount (sum ([droplets_size225] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize375ReturnCount (sum ([droplets_size375] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.
  let dropletSize750ReturnCount (sum ([droplets_size750] of patches with [returnVent = TRUE])) * ventil_movementRate  ; determine the number of droplets of this size that will be removed from returnVents.

  ask patches with [supplyVent = TRUE][ ; add droplets to supply vents from remove vents

    set additionalDroplets_size3  additionalDroplets_size3 + ((dropletSize3ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size6  additionalDroplets_size6 + ((dropletSize6ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size12 additionalDroplets_size12 + ((dropletSize12ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size20 additionalDroplets_size20 + ((dropletSize20ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size28 additionalDroplets_size28 + ((dropletSize28ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size36 additionalDroplets_size36 + ((dropletSize36ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size45 additionalDroplets_size45 + ((dropletSize45ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size62.5 additionalDroplets_size62.5 + ((dropletSize62.5ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size87.5 additionalDroplets_size87.5 + ((dropletSize87.5ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size112.5 additionalDroplets_size112.5 + ((dropletSize112.5ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size137.5 additionalDroplets_size137.5 + ((dropletSize137.5ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size175 additionalDroplets_size175 + ((dropletSize175ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size225 additionalDroplets_size225 + ((dropletSize225ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size375 additionalDroplets_size375 + ((dropletSize375ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.
    set additionalDroplets_size750 additionalDroplets_size750 + ((dropletSize750ReturnCount * (1 - ventil_removalRate)) / numSupplyVents)   ; The number of droplets added to supplyVents from returnVents is equal to the previously-summed number multiplied by the filtration rate divided by the number of supply vents, because the droplets are evenly distributed among supply vents.

  ]

  ask airArrows [ ; here we add droplets to patches due to air movement

    let dropletSize3Count [droplets_size3] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize6Count [droplets_size6] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize12Count [droplets_size12] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize20Count  [droplets_size20] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize28Count  [droplets_size28] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize36Count  [droplets_size36] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize45Count  [droplets_size45] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize62.5Count  [droplets_size62.5] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize87.5Count  [droplets_size87.5] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize112.5Count  [droplets_size112.5] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize137.5Count  [droplets_size137.5] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize175Count  [droplets_size175] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize225Count  [droplets_size225] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize375Count  [droplets_size375] of patch-here  ; determine the number of droplets of this size on the current patch.
  let dropletSize750Count  [droplets_size750] of patch-here  ; determine the number of droplets of this size on the current patch.

    ask patch-ahead 1 [ ;add droplets to patch 1 meter ahead at the pre-determined rate

    set additionalDroplets_size3  additionalDroplets_size3 + ( dropletSize3Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size6  additionalDroplets_size6 + ( dropletSize6Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size12 additionalDroplets_size12 + (dropletSize12Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size20 additionalDroplets_size20 + (dropletSize20Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size28 additionalDroplets_size28 + ( dropletSize28Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size36 additionalDroplets_size36 + ( dropletSize36Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size45 additionalDroplets_size45 + ( dropletSize45Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size62.5 additionalDroplets_size62.5 + ( dropletSize62.5Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size87.5 additionalDroplets_size87.5 + ( dropletSize87.5Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size112.5 additionalDroplets_size112.5 + ( dropletSize112.5Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size137.5 additionalDroplets_size137.5 + ( dropletSize137.5Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size175 additionalDroplets_size175 + ( dropletSize175Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size225 additionalDroplets_size225 + ( dropletSize225Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size375 additionalDroplets_size375 + ( dropletSize375Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.
    set additionalDroplets_size750 additionalDroplets_size750 + ( dropletSize750Count * ventil_movementRate)   ; The number of droplets added to forward patches is equal to the number of droplets on the current patch multiplied by the airflow rate.

    ]
  ]
  ;Now that we've denoted how many droplets of each size to add to patches, we can remove them. Note that this had to be done with a separate ask because if we add and remove in the same step, droplet additions to patches assessed later in the cue will be erroneous.

  ask patches [ ; remove droplets from all patches

    set droplets_size3 droplets_size3 - (droplets_size3 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size6 droplets_size6 - (droplets_size6 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size12 droplets_size12 - (droplets_size12 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size20 droplets_size20 - (droplets_size20 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size28 droplets_size28 - (droplets_size28 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size36 droplets_size36 - (droplets_size36 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size45 droplets_size45 - (droplets_size45 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size62.5 droplets_size62.5 - (droplets_size62.5 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size87.5 droplets_size87.5 - (droplets_size87.5 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size112.5 droplets_size112.5 - (droplets_size112.5 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size137.5 droplets_size137.5 - (droplets_size137.5 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size175 droplets_size175 - (droplets_size175 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size225 droplets_size225 - (droplets_size225 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size375 droplets_size375 - (droplets_size375 * ventil_movementRate) ; remove the appropriate number of droplets of this size
    set droplets_size750 droplets_size750 - (droplets_size750 * ventil_movementRate) ; remove the appropriate number of droplets of this size

  ]

  ;Now we officially add droplets, and update total droplet counts

    ask patches [ ; update droplet counts

    set droplets_size3 droplets_size3 + additionalDroplets_size3 ; add the appropriate number of droplets of this size
    set droplets_size6 droplets_size6 + additionalDroplets_size6 ; add the appropriate number of droplets of this size
    set droplets_size12 droplets_size12 + additionalDroplets_size12 ; add the appropriate number of droplets of this size
    set droplets_size20 droplets_size20 + additionalDroplets_size20 ; add the appropriate number of droplets of this size
    set droplets_size28 droplets_size28 + additionalDroplets_size28 ; add the appropriate number of droplets of this size
    set droplets_size36 droplets_size36 + additionalDroplets_size36 ; add the appropriate number of droplets of this size
    set droplets_size45 droplets_size45 + additionalDroplets_size45 ; add the appropriate number of droplets of this size
    set droplets_size62.5 droplets_size62.5 + additionalDroplets_size62.5 ; add the appropriate number of droplets of this size
    set droplets_size87.5 droplets_size87.5 + additionalDroplets_size87.5 ; add the appropriate number of droplets of this size
    set droplets_size112.5 droplets_size112.5 + additionalDroplets_size112.5 ; add the appropriate number of droplets of this size
    set droplets_size137.5 droplets_size137.5 + additionalDroplets_size137.5 ; add the appropriate number of droplets of this size
    set droplets_size175 droplets_size175 + additionalDroplets_size175 ; add the appropriate number of droplets of this size
    set droplets_size225 droplets_size225 + additionalDroplets_size225 ; add the appropriate number of droplets of this size
    set droplets_size375 droplets_size375 + additionalDroplets_size375 ; add the appropriate number of droplets of this size
    set droplets_size750 droplets_size750 + additionalDroplets_size750 ; add the appropriate number of droplets of this size

    set additionalDroplets_size3  0  ; reset the additional droplet counter
    set additionalDroplets_size6   0  ; reset the additional droplet counter
    set additionalDroplets_size12  0  ; reset the additional droplet counter
    set additionalDroplets_size20  0  ; reset the additional droplet counter
    set additionalDroplets_size28  0  ; reset the additional droplet counter
    set additionalDroplets_size36  0  ; reset the additional droplet counter
    set additionalDroplets_size45  0  ; reset the additional droplet counter
    set additionalDroplets_size62.5  0  ; reset the additional droplet counter
    set additionalDroplets_size87.5  0  ; reset the additional droplet counter
    set additionalDroplets_size112.5  0  ; reset the additional droplet counter
    set additionalDroplets_size137.5  0  ; reset the additional droplet counter
    set additionalDroplets_size175  0  ; reset the additional droplet counter
    set additionalDroplets_size225  0  ; reset the additional droplet counter
    set additionalDroplets_size375  0  ; reset the additional droplet counter
    set additionalDroplets_size750 0  ; reset the additional droplet counter

    set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count
    set virionCount (((1.40947E-14 * 1000 * virionsPerML) * droplets_size3) + ((1.12758E-13 * 1000 * virionsPerML) * droplets_size6) + ((9.02064E-13 * 1000 * virionsPerML) * droplets_size12) + ((4.17622E-12 * 1000 * virionsPerML) * droplets_size20) + ((1.14595E-11 * 1000 * virionsPerML) * droplets_size28) + ((2.43557E-11 * 1000 * virionsPerML) * droplets_size36) + ((4.75698E-11 * 1000 * virionsPerML) * droplets_size45) + ((1.27448E-10 * 1000 * virionsPerML) * droplets_size62.5) + ((3.49718E-10 * 1000 * virionsPerML) * droplets_size87.5) + ((7.43277E-10 * 1000 * virionsPerML) * droplets_size112.5) + ((1.35707E-09 * 1000 * virionsPerML) * droplets_size137.5) + ((2.79774E-09 * 1000 * virionsPerML) * droplets_size175) + ((5.94622E-09 * 1000 * virionsPerML) * droplets_size225) + ((2.75288E-08 * 1000 * virionsPerML) * droplets_size375) + ((2.2023E-07 * 1000 * virionsPerML) * droplets_size750)) ; Multiply estimated virion counts per droplet by droplet counts (Note that the first term is in each equation is the mass or volume of each droplet size in kg or L (it's interchangeable because the weight of 1 liter of pure water at temperature 4 °C = 1) . We multiply this by 1000 to get the volume of each droplet size in mL.
    set transmissionRisk (virionCount * virionRisk) ;update risk of infection in patches.

  ]
  ]

end

to removeDroplets ; Here we remove droplets due to: inhalation by people, gravitational settling, and decay.
  ; Note that when removing droplets from patches, we make several assumptions.
  ; 1.) "droplet removal" is equivalent to virion removal.
  ; 2.) All droplet sizes are post-evaporation, so there are no evaporation rates included in removal equations.
  ; 3.) The decay rate is independent of droplet size.
  ; 4.) The inhalation rate is constant for all people.
  ; 5.) Gravitational fallout rate (per minute) is based on droplet diameter (in micrometers) and is calculated using the equations described by Anchordoqui & (Chudnovsky 2020).
  ;     Anchordoqui, L.A., & Chudnovsky, E.M. 2020. A physicist view of COVID-19 airborne infection through convective airflow in indoor spaces. Preprint available at https://arxiv.org/pdf/2003.13689.pdf.
  ;     We then convert the time to ground-level (from expectorate-height) estimate for droplets of different diameters from seconds to rate per minute.

  ask patches with [totalDroplets > 0][ ; no need to ask ALL patches, just the ones with droplets present

    let inhaledPercent (vol_B / (1 * 1 * expectorate-height)) * (count people-here) ; this is the fraction of air in the patch inhaled by individuals each minute. We assume that droplets are well mixed and evenly distributed within patches, that people only have access to droplet-height cubic meters of air, and that the breathing rate is constant for all individuals.
    let size3_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam3) / 60))  + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size3_change > 1) [set size3_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size6_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam6) / 60)) + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size6_change > 1) [set size6_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size12_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam12) / 60)) + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size12_change > 1) [set size12_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size20_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam20) / 60)) + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size20_change > 1) [set size20_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size28_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam28) / 60)) + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size28_change > 1) [set size28_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size36_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam36) / 60)) + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size36_change > 1) [set size36_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size45_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam45) / 60)) + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size45_change > 1) [set size45_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size62.5_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam62.5) / 60)) + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size62.5_change > 1) [set size62.5_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size87.5_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam87.5) / 60))  + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size87.5_change > 1) [set size87.5_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size112.5_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam112.5) / 60))  + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size112.5_change > 1) [set size112.5_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size137.5_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam137.5) / 60))  + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size137.5_change > 1) [set size137.5_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size175_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam175) / 60))  + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size175_change > 1) [set size175_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size225_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam225) / 60))  + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size225_change > 1) [set size225_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size375_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam375) / 60))  + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size375_change > 1) [set size375_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.
    let size750_change (inhaledPercent + (1 / ((expectorate-height / Vt_diam750) / 60))  + dropletDecay) ; These rates will be added together to determine the proportion of droplets removed each tick.
    if (size750_change > 1) [set size750_change 1] ; If the equation says the removed proportion is greater than 100%, we cap it at 1.

    ;Update droplet totals

    set droplets_size3 (droplets_size3 - (size3_change * droplets_size3)) ;remove the determined proportion of size 3 droplets.
    set droplets_size6 (droplets_size6 - (size6_change * droplets_size6)) ;remove the determined proportion of size 6 droplets.
    set droplets_size12 (droplets_size12 - (size12_change * droplets_size12)) ;remove the determined proportion of size 12 droplets.
    set droplets_size20 (droplets_size20 - (size20_change * droplets_size20)) ;remove the determined proportion of size 20 droplets.
    set droplets_size28 (droplets_size28 - (size28_change * droplets_size28)) ;remove the determined proportion of size 28 droplets.
    set droplets_size36 (droplets_size36 - (size36_change * droplets_size36)) ;remove the determined proportion of size 36 droplets.
    set droplets_size45 (droplets_size45 - (size45_change * droplets_size45)) ;remove the determined proportion of size 45 droplets.
    set droplets_size62.5 (droplets_size62.5 - (size62.5_change * droplets_size62.5)) ;remove the determined proportion of size 62.5 droplets.
    set droplets_size87.5 (droplets_size87.5 - (size87.5_change * droplets_size87.5)) ;remove the determined proportion of size 87.5 droplets.
    set droplets_size112.5 (droplets_size112.5 - (size112.5_change * droplets_size112.5)) ;remove the determined proportion of size 112.5 droplets.
    set droplets_size137.5 (droplets_size137.5 - (size137.5_change * droplets_size137.5)) ;remove the determined proportion of size 137.5 droplets.
    set droplets_size175 (droplets_size175 - (size175_change * droplets_size175)) ;remove the determined proportion of size 175 droplets.
    set droplets_size225 (droplets_size225 - (size225_change * droplets_size225)) ;remove the determined proportion of size 225 droplets.
    set droplets_size375(droplets_size375 - (size375_change * droplets_size375)) ;remove the determined proportion of size 375 droplets.
    set droplets_size750 (droplets_size750 - (size750_change * droplets_size750)) ;remove the determined proportion of size 750 droplets.

    set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count
    set virionCount (((1.40947E-14 * 1000 * virionsPerML) * droplets_size3) + ((1.12758E-13 * 1000 * virionsPerML) * droplets_size6) + ((9.02064E-13 * 1000 * virionsPerML) * droplets_size12) + ((4.17622E-12 * 1000 * virionsPerML) * droplets_size20) + ((1.14595E-11 * 1000 * virionsPerML) * droplets_size28) + ((2.43557E-11 * 1000 * virionsPerML) * droplets_size36) + ((4.75698E-11 * 1000 * virionsPerML) * droplets_size45) + ((1.27448E-10 * 1000 * virionsPerML) * droplets_size62.5) + ((3.49718E-10 * 1000 * virionsPerML) * droplets_size87.5) + ((7.43277E-10 * 1000 * virionsPerML) * droplets_size112.5) + ((1.35707E-09 * 1000 * virionsPerML) * droplets_size137.5) + ((2.79774E-09 * 1000 * virionsPerML) * droplets_size175) + ((5.94622E-09 * 1000 * virionsPerML) * droplets_size225) + ((2.75288E-08 * 1000 * virionsPerML) * droplets_size375) + ((2.2023E-07 * 1000 * virionsPerML) * droplets_size750)) ; Multiply estimated virion counts per droplet by droplet counts (Note that the first term is in each equation is the mass or volume of each droplet size in kg or L (it's interchangeable because the weight of 1 liter of pure water at temperature 4 °C = 1) . We multiply this by 1000 to get the volume of each droplet size in mL.
    set transmissionRisk (virionCount * virionRisk) ;update risk of infection in patches.
  ]

end

to testQuanta ;determine the number of people that will be infected if the number of droplets was distributed over the entire world

  ;First we contaminate patches

    let will_cough random-float 1 ; pull a number random number between 0 and 1 to determine if people will cough, or if the expectoration event will be parameterized to reflect speaking.

    ifelse(will_cough <= cough_frequency)[ ; if will_cough is less than or equal to cough_frequency, the resulting cone of exposure will be produced based on cough parameters. If not, the cone of exposure will be based on speaking parameters.

          let dropletNum exp (random-normal M_cough_droplets S_cough_droplets) ; pull the number of droplets that will be added to patch counts from a lognormal distribution following the procedure described by (Railsback & Grimm 2011)

        ask patches [ ; first we equally distribute droplets across all patches
          set droplets_size3 (droplets_size3 + (dropletNum * 0.091863517))  ;Add the number of droplets multiplied by the probability that droplets will be 2-4 (i.e., mean of 3) micrometers.
          set droplets_size6 (droplets_size6 + (dropletNum * 0.461942257))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + (dropletNum * 0.170603675))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + (dropletNum * 0.073490814))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.
          set droplets_size28 (droplets_size28 + (dropletNum * 0.036745407))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.

      ;at a start height of 1.7m, all droplets in size clases with diameters larger than 28 micrometers fall out in less than 1 minute

          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count
          set virionCount (((1.40947E-14 * 1000 * virionsPerML) * droplets_size3) + ((1.12758E-13 * 1000 * virionsPerML) * droplets_size6) + ((9.02064E-13 * 1000 * virionsPerML) * droplets_size12) + ((4.17622E-12 * 1000 * virionsPerML) * droplets_size20) + ((1.14595E-11 * 1000 * virionsPerML) * droplets_size28) + ((2.43557E-11 * 1000 * virionsPerML) * droplets_size36) + ((4.75698E-11 * 1000 * virionsPerML) * droplets_size45) + ((1.27448E-10 * 1000 * virionsPerML) * droplets_size62.5) + ((3.49718E-10 * 1000 * virionsPerML) * droplets_size87.5) + ((7.43277E-10 * 1000 * virionsPerML) * droplets_size112.5) + ((1.35707E-09 * 1000 * virionsPerML) * droplets_size137.5) + ((2.79774E-09 * 1000 * virionsPerML) * droplets_size175) + ((5.94622E-09 * 1000 * virionsPerML) * droplets_size225) + ((2.75288E-08 * 1000 * virionsPerML) * droplets_size375) + ((2.2023E-07 * 1000 * virionsPerML) * droplets_size750)) ; Multiply estimated virion counts per droplet by droplet counts (Note that the first term is in each equation is the mass or volume of each droplet size in kg or L (it's interchangeable because the weight of 1 liter of pure water at temperature 4 °C = 1) . We multiply this by 1000 to get the volume of each droplet size in mL.
          set transmissionRisk (virionCount * virionRisk) ;update risk of infection in patches.

    ]
   ]

  [ ; if we're not drawing from a coughing distribution (i.e., we're drawsing from a speaking one.
    let dropletNum exp (random-normal M_speak_droplets S_speak_droplets) ; pull the number of droplets that will be added to patch counts from a lognormal distribution following the procedure described by (Railsback & Grimm 2011)

  ask patches [ ; first we equally distribute droplets across all patches

          set droplets_size3 (droplets_size3 + ((dropletNum * 0.091863517) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 2-4 (i.e., mean of 3) micrometers.
          set droplets_size6 (droplets_size6 + ((dropletNum * 0.461942257) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 4-8 (i.e., mean of 6) micrometers.
          set droplets_size12 (droplets_size12 + ((dropletNum * 0.170603675) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 8-16 (i.e., mean of 12) micrometers.
          set droplets_size20 (droplets_size20 + ((dropletNum * 0.073490814) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 16-24 (i.e., mean of 20) micrometers.

        ;at a start height of 1.7m, all droplets in size clases with diameters larger than 20 micrometers fall out in less than 2 minutes

          set droplets_size28 (droplets_size28 + ((dropletNum * 0.036745407) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 24-32 (i.e., mean of 28) micrometers.

    ;at a start height of 1.7m, all droplets in size clases with diameters larger than 28 micrometers fall out in less than 1 minute

;          set droplets_size36 (droplets_size36 + ((dropletNum * 0.015748031) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 32-40 (i.e., mean of 36) micrometers.
;          set droplets_size45 (droplets_size45 + ((dropletNum * 0.005249344) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 40-50 (i.e., mean of 45) micrometers.
;          set droplets_size62.5 (droplets_size62.5 + ((dropletNum * 0.023622047) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 50-75 (i.e., mean of 62.5) micrometers.
;          set droplets_size87.5 (droplets_size87.5 + ((dropletNum * 0.01312336) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 75-100 (i.e., mean of 87.5) micrometers.
;          set droplets_size112.5 (droplets_size112.5 + ((dropletNum * 0.026246719) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 100-125 (i.e., mean of 112.5) micrometers.
;          set droplets_size137.5 (droplets_size137.5 + ((dropletNum * 0.018372703) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 125-150 (i.e., mean of 137.5) micrometers.
;          set droplets_size175 (droplets_size175 + ((dropletNum * 0.015748031) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 150-200 (i.e., mean of 175) micrometers.
;          set droplets_size225 (droplets_size225 + ((dropletNum * 0.01312336) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 200-250 (i.e., mean of 225) micrometers.
;          set droplets_size375 (droplets_size375 + ((dropletNum * 0.023622047) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 250-500 (i.e., mean of 375) micrometers.
;          set droplets_size750 (droplets_size750 + ((dropletNum * 0.010498688) / (count patches)))  ;Add the number of droplets multiplied by the probability that droplets will be 500-1000 (i.e., mean of 750) micrometers.
          set totalDroplets (droplets_size3 + droplets_size6 + droplets_size12 + droplets_size20 + droplets_size28 + droplets_size36 + droplets_size45 + droplets_size62.5 + droplets_size87.5 + droplets_size112.5 + droplets_size137.5 + droplets_size175 + droplets_size225 + droplets_size375 + droplets_size750) ; update the total droplet count
          set virionCount (((1.40947E-14 * 1000 * virionsPerML) * droplets_size3) + ((1.12758E-13 * 1000 * virionsPerML) * droplets_size6) + ((9.02064E-13 * 1000 * virionsPerML) * droplets_size12) + ((4.17622E-12 * 1000 * virionsPerML) * droplets_size20) + ((1.14595E-11 * 1000 * virionsPerML) * droplets_size28) + ((2.43557E-11 * 1000 * virionsPerML) * droplets_size36) + ((4.75698E-11 * 1000 * virionsPerML) * droplets_size45) + ((1.27448E-10 * 1000 * virionsPerML) * droplets_size62.5) + ((3.49718E-10 * 1000 * virionsPerML) * droplets_size87.5) + ((7.43277E-10 * 1000 * virionsPerML) * droplets_size112.5) + ((1.35707E-09 * 1000 * virionsPerML) * droplets_size137.5) + ((2.79774E-09 * 1000 * virionsPerML) * droplets_size175) + ((5.94622E-09 * 1000 * virionsPerML) * droplets_size225) + ((2.75288E-08 * 1000 * virionsPerML) * droplets_size375) + ((2.2023E-07 * 1000 * virionsPerML) * droplets_size750)) ; Multiply estimated virion counts per droplet by droplet counts (Note that the first term is in each equation is the mass or volume of each droplet size in kg or L (it's interchangeable because the weight of 1 liter of pure water at temperature 4 °C = 1) . We multiply this by 1000 to get the volume of each droplet size in mL.
          set transmissionRisk (virionCount * virionRisk) ;update risk of infection in patches.

  ]
  ]

        ;Now that we have contaminated patches, we can infect others. Infection probability is determined by the cumulative sum of virions in droplets in the patch multiplied by virionRisk (i.e., the risk of infection given exposure to a single virion).
      ; For SARS-CoV2, we assume there are 7 X 10^6 functional virions per milliliter of droplet fluid (Wolfel et al. 2020). We calculate the mass (in kg) of each droplet using the equation described by Anchordoqui & Chudnovski (2020), and convert kg to L, then to mL for virion estimation.
      ; Additionally, we assume that droplets are evenly distributed within patches, and that the air within which droplets are distubited is contained within 1 m X 1 m X expectorate-height m cubic area (i.e., patches represent 1 m X 1 m X expectorate-height m cubic space, and individuals within these patches must breathe in air from this space).

  color-patches

   ask patches [ ; call on infectious patches.

      ask people-here [ ;ask people in infectious patches

        set timesteps_exposed (timesteps_exposed + 1) ; add 1 to the timesteps_exposed variable (i.e., counts the number of times individuals were on infectious patches)

        if(infected? = FALSE AND infectious? = FALSE)[ ; only uninfected people can become infected

           let will_infect? random-float 1 ; pull a random number between zero and 1
           if(will_infect? <= (exposure-risk * ([transmissionRisk] of patch-here)  * (vol_B / (1 * 1 * expectorate-height)) ))[ ; if will_infect? does not exceed the risk of infection (i.e., the probability that a person is exposed to virions in the patch AND that they breathe it in AND that inhilation will lead to infection), the person is infected.

                set color 25 ; infected people will be this shade of orange
                set infected? TRUE ; update infection status
                set totalInfected (totalInfected + 1) ; add to the totalInfected count

              if(totalInfected = 1)[ ; when the totalInfected count first becomes 1, we record the time at which it happens

               set firstInfectTime ticks ; record the time of the first infection event.

              ]


            ]
            ]
    ]
    ]

  show totalInfected / (num_completelySusceptible + num_reducedSusceptible)

  stop

end
@#$#@#$#@
GRAPHICS-WINDOW
1338
292
1706
661
-1
-1
40.0
1
10
1
1
1
0
0
0
1
0
8
0
8
0
0
1
ticks
30.0

INPUTBOX
26
196
99
256
grid-height
9.0
1
0
Number

INPUTBOX
99
196
173
256
grid-width
9.0
1
0
Number

BUTTON
26
58
92
91
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
92
58
158
91
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
158
58
224
91
NIL
clear-all
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
24
427
98
487
n
135.0
1
0
Number

INPUTBOX
248
541
338
601
social-distance
0.0
1
0
Number

MONITOR
177
204
250
249
n patches
count patches
17
1
11

INPUTBOX
25
287
91
347
cohort-dur
60.0
1
0
Number

INPUTBOX
248
481
338
541
maskRisk-mod
0.1
1
0
Number

INPUTBOX
655
221
784
281
cough_spread-dist.mean
5.0
1
0
Number

INPUTBOX
784
221
901
281
cough_spread-dist.sd
0.256
1
0
Number

PLOT
1331
57
1498
213
Total Infected
time (mins)
infected people
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot totalInfected"

TEXTBOX
1812
11
2454
206
MODEL OUTPUTS: (these are the global variables that can be retrieved for analysis following simulations)\n\nfirstInfectTime - Time (i.e., tick) at which the first successful transmission event occurs. \n\ntotalInfected    -    Total number of infected people over the course of the simulation (i.e., running sum).\n\ntotalInfected.list    -    Total number of infected people at each tick.\n\navg.dist - Average interpersonal distance at a given tick.\n\navg.exposures - Average count of instances when susceptible agents were within a cone of exposure.\n\navg.exposureTime - Average number of timesteps when susceptible agents were within ≥ 1 cone of exposure. \n\n\n
10
0.0
1

INPUTBOX
548
221
656
281
cough_airflow-angle
35.0
1
0
Number

MONITOR
1330
12
1562
57
Number of symptomatic individuals
num_symptomatic
17
1
11

INPUTBOX
548
280
662
340
speak_airflow-angle
63.5
1
0
Number

INPUTBOX
661
280
797
340
speak_spread-dist.mean
0.55
1
0
Number

INPUTBOX
796
280
915
340
speak_spread-dist.sd
0.068
1
0
Number

TEXTBOX
1513
63
1746
258
COLOR LEGEND\n\nPeople:\n\nRed - symptomatic\nPink - asymptomatic\nGreen - susceptible\nOrange - infected\n\nPatches:\n\nCyan - uncontaminated\nMagenta - droplet contaminated. Darkness scales with level of infectiousness. Incredibly contaminated patches will appear black\n
10
0.0
1

INPUTBOX
248
422
339
482
mod-proportion
0.0
1
0
Number

INPUTBOX
750
161
845
221
cough_frequency
0.19
1
0
Number

INPUTBOX
24
487
98
547
n_infectious
1.0
1
0
Number

INPUTBOX
24
546
98
606
symp-pr
1.0
1
0
Number

MONITOR
1561
12
1800
57
Number of asymptomatic individuals
num_asymptomatic
17
1
11

CHOOSER
550
583
671
628
mod_group
mod_group
"sus" "inf" "sus_inf"
2

INPUTBOX
548
340
680
400
speak_dropletNum.mean
1000.0
1
0
Number

INPUTBOX
679
340
799
400
speak_dropletNum.sd
0.0
1
0
Number

INPUTBOX
797
340
933
400
cough_dropletNum.mean
1000.0
1
0
Number

INPUTBOX
933
340
1054
400
cough_dropletNum.sd
0.0
1
0
Number

INPUTBOX
241
1077
318
1137
virionRisk
0.0624
1
0
Number

INPUTBOX
230
600
338
660
expectorate-height
1.7
1
0
Number

INPUTBOX
24
605
98
665
vol_B
0.023
1
0
Number

INPUTBOX
1016
160
1104
220
dropletDecay
0.9
1
0
Number

SWITCH
29
682
183
715
face-northward
face-northward
1
1
-1000

INPUTBOX
550
493
662
553
personPerPatch-cap
2.0
1
0
Number

INPUTBOX
240
284
316
344
num-cohorts
1.0
1
0
Number

SWITCH
26
359
171
392
rearrange-cohort
rearrange-cohort
1
1
-1000

INPUTBOX
29
1076
104
1136
virionsPerML
2.35E9
1
0
Number

SWITCH
29
745
150
778
ventilation
ventilation
0
1
-1000

INPUTBOX
28
777
127
837
ventilSupplyWall
east
1
0
String

INPUTBOX
126
777
225
837
ventilReturnWall
west
1
0
String

INPUTBOX
28
836
126
896
numSupplyVents
3.0
1
0
Number

INPUTBOX
125
836
225
896
numReturnVents
3.0
1
0
Number

INPUTBOX
28
926
153
986
ventil_movementRate
0.5
1
0
Number

INPUTBOX
28
986
141
1046
ventil_removalRate
1.0
1
0
Number

SWITCH
28
895
174
928
equallySpaceVents
equallySpaceVents
0
1
-1000

INPUTBOX
548
162
631
222
diffusionRate
0.0015
1
0
Number

SWITCH
556
728
682
761
evalQuanta
evalQuanta
1
1
-1000

MONITOR
1361
76
1418
121
num
totalInfected
0
1
11

SWITCH
550
552
671
585
showArrows
showArrows
0
1
-1000

CHOOSER
548
399
714
444
speak_dropletSizeDistr
speak_dropletSizeDistr
"chao" "meanlog.1" "meanlog.2" "meanlog.3" "meanlog.4" "meanlog.5"
0

CHOOSER
714
399
881
444
cough_dropletSizeDistr
cough_dropletSizeDistr
"chao" "meanlog.1" "meanlog.2" "meanlog.3" "meanlog.4" "meanlog.5"
0

TEXTBOX
27
12
92
35
Controls 
14
0.0
1

TEXTBOX
26
104
253
138
Scenario-specific model inputs
14
0.0
1

TEXTBOX
28
35
480
63
Initialize, run, or clear the simulation.
11
0.0
1

TEXTBOX
25
126
446
168
These parameters describe agent/environmental behavior in specific circumstances. Users are likely to change them from their current values. Input categories are shown in orange.
11
0.0
1

TEXTBOX
254
199
404
254
Control height and width of the simulated world. The \"n patches\" reporter shows the number of patches (i.e., square meters) in the simulated area.
9
0.0
1

TEXTBOX
104
439
237
472
Total number of individuals (susceptible + infectious) in the simulation.
9
0.0
1

TEXTBOX
104
507
254
535
Number of infectious individuals.
9
0.0
1

TEXTBOX
105
559
232
592
The probability (0 - 1) that infectious individuals will be symptomatic.
9
0.0
1

TEXTBOX
25
172
175
190
World size
13
24.0
1

TEXTBOX
26
402
223
434
Individual characteristics
13
24.0
1

TEXTBOX
25
263
175
281
Simulation length
13
24.0
1

TEXTBOX
95
299
217
332
Number of ticks each cohort of individuals spends in the simulation.
9
0.0
1

TEXTBOX
322
299
472
321
Number of cohorts to be simulated.
9
0.0
1

TEXTBOX
179
352
450
402
Controls whether or not cohorts are replaced after \"cohort-dur\" ticks up to \"num-cohorts\" times, or if the locations of individuals in the original cohort just have their locations rearranged \"num-cohorts\" times.
9
0.0
1

TEXTBOX
346
435
496
468
The probability (0 - 1) that individuals will be wearing masks.
9
0.0
1

TEXTBOX
345
556
468
596
The distance (m) that individuals attempt to maintain from others.
9
0.0
1

TEXTBOX
190
681
340
714
Controls if individuals all face north (upwards), or if they face a random direction.
9
0.0
1

TEXTBOX
30
722
180
740
Ventilation
13
24.0
1

TEXTBOX
155
750
320
783
Dictates whether or not to simulate forced-air ventilation effects.
9
0.0
1

TEXTBOX
345
501
472
534
The proportion of droplets (0 - 1) exhaled/inhaled by mask wearers.
9
0.0
1

TEXTBOX
678
551
881
588
Control if users want to see arrows denoting the direction of airflow when simulating ventilation effects.
9
0.0
1

TEXTBOX
233
795
383
817
Controls the locations of any supply and return vents.
9
0.0
1

TEXTBOX
234
859
384
881
Controls the number of any supply and return vents.
9
0.0
1

TEXTBOX
182
896
375
940
Control whether vent locations are equally spaced along the walls, or if their locations are randomly decided.
9
0.0
1

TEXTBOX
162
943
364
998
Rate of droplet movement (% droplets / tick) between patches when ventilation effects are simulated. Values range from 0 - 1.
9
0.0
1

TEXTBOX
148
993
298
1037
Percent (0 - 1) of droplets removed from the simulation due to filtration at return-vent patches every tick. 
9
0.0
1

TEXTBOX
29
1054
179
1072
Pathogen
13
24.0
1

TEXTBOX
110
1095
232
1117
Number of virions present per mL of droplet fluid. 
9
0.0
1

TEXTBOX
326
1086
476
1142
The probability that inhalation of a virion will result in infection for a susceptible individual.
9
0.0
1

TEXTBOX
548
18
698
36
Stable inputs
14
0.0
1

TEXTBOX
548
39
835
129
These parameters generally reflect real-world droplet dynamics and other parameters that are less likely to be changed, regardless of the simulated scenario. See Farthing et al. 2021 - the paper wherein we introduced this model, for a detailed description of value determination. Categories are shown in red.
11
0.0
1

TEXTBOX
637
173
744
210
Rate (%/tick) at which droplets move to neighboring patches.
9
0.0
1

TEXTBOX
670
504
820
537
Maximum number of people allowed in a single patch at any given time.
9
0.0
1

TEXTBOX
678
592
840
625
Dictates whether infectious people only, susceptible people only, or both wear masks. 
9
0.0
1

TEXTBOX
106
623
201
655
Inhalation rate (cubic meters / tick)
9
0.0
1

TEXTBOX
345
613
495
657
Height of individuals (m) within the simulation (represents height from which expectorate falls).
9
0.0
1

TEXTBOX
548
137
698
158
Droplet dynamics
13
13.0
1

TEXTBOX
550
464
700
482
Other
13
13.0
1

TEXTBOX
852
172
1002
205
Probabilty (0 - 1) that symptomatic infectious individuals will cough each tick. 
9
0.0
1

TEXTBOX
907
227
1057
271
Angle (degree), mean and standard-deviation distance (m) of droplet cone spread when coughing.  
9
0.0
1

TEXTBOX
921
281
1071
336
Angle (degree), mean and standard-deviation distance (m) of droplet cone spread when NOT coughing (parameterized as \"speaking\" events).  
9
0.0
1

TEXTBOX
1061
345
1211
389
Number of droplets expelled during coughing and non-coughing (i.e., \"speaking\") events.
9
0.0
1

TEXTBOX
889
402
1213
531
Defines the probability distribution for droplet sizes during coughing and non-coughing (i.e., \"speaking\") events. If \"chao,\" the distribution is defined using the findings of Chao et al. (2009). Otherwise, droplet sizes are drawn from a lognormal distribution with a log standard deviation of 3 micrometers, and a log mean defined here. \n\nChao, C.Y.H., Wan, M.P., Morawska, L., Johnson, G.R., Ritovski, Z.D., …, & Katoshevski, D. (2009). Characterization of expiration air jets and droplet size distributions immediately at the mouth opening. Aerosol Science 40(2009):122 – 133. https://doi.org/10.1016/j.jaerosci.2008.10.003.\n
9
0.0
1

TEXTBOX
1110
176
1204
198
Droplet/virion decay rate (% / tick).
9
0.0
1

TEXTBOX
553
651
703
669
Special
14
0.0
1

TEXTBOX
556
674
834
722
These are special-case inputs that instruct our model to run sub-models that fundamentally changes the way the model works.
11
0.0
1

TEXTBOX
688
727
1259
780
If \"evalQuanta\" is turned on, the droplet fallout procedure will be carried out before the infection sub-model (so that only aerosolized droplets can trigger infection), and droplets are homogenously dispersed throughout the entirety of the simulated world immediately after expectoration. This can be used to find the combination of pathogen/droplet dynamics parameters that effectively recreate one quantum (i.e., the number of aerosolized infectious particles required to infect 1- 1/e % of a population) in a specific scenario. Currently evalQuanta is only coded to allow for an expectorate-height of 1.7 m.
9
0.0
1

SWITCH
556
799
682
832
choirSim
choirSim
1
1
-1000

TEXTBOX
689
798
1256
886
\"choirSim\" is used to simulate the SARS-CoV-2 superspreading event that took place dascribed by Hamner et al. (2020). If turned on, the individual characteristics, simulation length, and world size scenario-specific inputs are forcefully changed to reflect this event. Additionally, as the simulation runs, the cohort of individuals is rearranged 4 times, on the fourth time they return to their original positions. \n\nHamner, L., Dubbel, P., Capron, I., Ross, A., Jordan, A., Lee, J., …, & Leibrand, H. (2020). High SARS-CoV-2 attack rate following exposure at a choir practice – Skagit County, Washington, March 2020. Morbidity and Mortality Weekly Report 69(2020):606-610. http://dx.doi.org/10.15585/mmwr.mm6919e6.
9
0.0
1

@#$#@#$#@
## WHAT IS IT?

This is a spatially-explicit, stochastic agent-based model for simulating airborne and direct droplet-mediated respiratory pathogen transmission in indoor settings. Its purpose is to quantify the effect of increasing group density on the probability of respiratory pathogen transmission from infectious individuals to susceptible ones, given varied spatial dimension and risk-reduction behavior (e.g., mask use, social distancing, etc.) levels in indoor settings. The ability of this model to accurately simulate infection events is predicated on its ability to recreate four processes involved in transmission: 1.) Susceptible individuals become infected through inhalation of virions contained within infectious droplets of varying sizes. 2.) Infectious agents expel infectious droplets of varying sizes, and droplets’ movement, fallout, and virion-carriage rates vary with droplet size. 3.) Symptomatic infectious agents are likely to infect more susceptible individuals than asymptomatic ones, as coughing expels infectious droplets farther than does breathing or speaking alone. 4.) Susceptible individuals’ probability of infection can be lessened if individuals employ extra measures to avoid transmission (e.g., wearing face masks).

## HOW IT WORKS

For a detailed description of how the model works, see our Overview, Design concepts, Details (ODD) model description available at: and on the Lanzaslab github page (https://github.com/lanzaslab/droplet-ABM).

## HOW TO USE IT

All global variables are described in detail on the CODE tab. Users may change model input values on the interface tab to reflect activity- or case-specific scenarios, then run the simulation to assess how respiratory-pathogen infections may propagate in indoor setttings. Model output variables are described in a note on the interface tab. In Farthing et al. 2021 we use this model to estimate effects of proposed COVID-19 intervention strategies for indoor environments (i.e., increased airflow, limiting contact durations, wearing masks, and increased interpersonal spacing), and investigate potential drivers of airborne superspreading events.

## EXTENDING THE MODEL

There are a number of areas that we plan to improve in future iterations of this model, but users may code them in themselves if they so choose (e.g., more complex movement capabilities, the ability to more-effectively recreate specific activities like dining in a restaurant, more-complex ventilation systems, etc.).


## CREDITS AND REFERENCES

When referencing this model in any presentations or publications, users should cite:
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
