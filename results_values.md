## Optimally doped YBCO, HPHT

|  Optimization and model    | Temperature range | β (en mJ/mol/K⁴) | Θ (en K) | γ (en mJ/mol/K²) | Remarques                                        |
|--------------------------|----------------------|------------------|----------|-------------------|--------------------------------------------------|
| Debye non linear ($Cs + \beta x + \gamma,  x=T²$)  | 0 - 10 K             | 0.486            | 373      | 1.55           | E seems under-estimated, values consistent with those found with Monte Carlo method, $\chi ² = 8.86$ |
| Debye linear  ($\beta x + \gamma,  x=T^2$) | 6 - 12 K             | 0.481            | 375      | 4.43         | Consistent with Monte Carlo, $\chi ² = 12.13$|
| Debye linear ($\beta x + \gamma,  x=T^2$) | 5 - 12 K             | 0.453            | 382      | 7.05       |  Consistent with Monte Carlo, $\chi ² = 50.59$  |
| Contribution en T⁵, non linéaire ($\beta x + \gamma + \alpha x² + Cs,  x=T^2$)  | 0 - 20 K             | 0.459            | 380      | 0.566             | Consistent with MC, $\chi ² = 63.78$ |
| Contribution en T⁵ polynomial ($\beta x + \gamma + \alpha x²$) | 6 - 20 K            | 0.450            | 388      | 1.01             |  Same results with MC, $\chi ² = 80$  |

Note : The linear fit performed between 5 and 12 K might not be completely free of the Schottky anomaly, since the region in which it can be neglected has been determined with a parameter E likely underestimated, so its maximum is shifted to the left. Morevover, we can see the $\chi ²$ is really lower for the 6-12 K fit.

## Overdoped YBCO, ISIS sample

| Optimization, model    | Temperature range | β (en mJ/mol/K⁴) | Θ (en K) | γ (en mJ/mol/K²) | Remarques                                        |
|--------------------------|----------------------|------------------|----------|-------------------|--------------------------------------------------|
| Debye non linear ($Cs + \beta x + \gamma,  x=T²$)  | 0 - 12 K             | 0.750            | 324      | 3.16              |same values with MC, $\chi ²$ = 47|
| Debye linear   ($\beta x + \gamma,  x=T^2$) | 5 - 12 K             | 0.772            | 327      | 6.95              | Matches Monte Carlo $\chi ²$ = 72    |
| Debye linear   ($\beta x + \gamma,  x=T^2$)| 5 - 10 K             | 0.676            | 334      | 9.00  |  $\chi ²$ = 168 |     
| Contribution T⁵ non linear ($\beta x + \gamma + \alpha x² + Cs,  x=T^2$)  | 0 - 20 K             | 0.738            | 324      | 2.31              | Same values with MC, $\chi ²$ = 202 |
| Contribution T⁵ polynomial ($\beta x + \gamma + \alpha x²$) | 5 - 20 K             | 0.711            | 328      | 5.73              |   Matches the values found with MC, $\chi ²$ = 296 |

## Overdoped YBCO, P240917 sample

| Optimization, model     | Temperature range | 𝛽 (en mJ/mol/K⁴) | Θ (en K) | γ (en mJ/mol/K²) | Remarques                                  |
| :----------------------- | :------------------- | :--------------- | :------- | :--------------- | :----------------------------------------- |
| Debye non linear       | 0 - 12 K             | 0.910            | 303      | 10.95            |   Matches the MC values, $\chi ² = 43$    |
| Debye, non linear     | 0 - 20 K             | 0.992            | 294      | 2.99             |Did not work for the others samples|
| Debye linear           | 5 - 12 K             | 0.885            | 305      | 14.26            |   Matches the MC values, $\chi ² = 57$  |
| Debye linear          | 5 - 10 K             | 0.807            | 314      | 18.22            |   Matches the MC values, $\chi ² = 10.26$    |
| Contribution T⁵ non linear    | 0 - 20 K             | 0.910            | 302      | 9.32 | Matches the MC values, $\chi ² = 97$   | 
| Contribution T⁵ polynomial      | 5 - 20 K             | 0.879            | 306      | 13.37         |     Matches the MC values, $\chi ² = 134$     |