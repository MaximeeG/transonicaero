# Q2 report figures

P01–P20: density, pressure, velocity, and Mach profiles.
C01–C20: density convergence histories.
G01–G08: paired profile/convergence comparisons for the CFL and grid studies.

Each CSV is read once. Re-running regenerates these figures.

Profiles are final saved states. Invalid profiles are replaced by explanatory panels.
Convergence is log10 of the unnormalized L2 density change between consecutive saved iterations.
Stopping status is reconstructed from rounded CSV values with tolerance 1e-8; it is not a solver-certified status.
No gamma or dissipation settings are inferred from filenames.

| Run | Settings | Outcome | Figures |
|---|---|---|---|
| 1 | MacCormack ; supersonic outlet ; CFL = 0.1 ; Nx = 201 | not converged at final saved iteration | [P01](P01_MACCORMACK_SUPER_CFL_10_NX_201.png), [C01](C01_MACCORMACK_SUPER_CFL_10_NX_201.png) |
| 2 | MacCormack ; supersonic outlet ; CFL = 0.5 ; Nx = 201 | near stopping tolerance (CSV precision limited) | [P02](P02_MACCORMACK_SUPER_CFL_50_NX_201.png), [C02](C02_MACCORMACK_SUPER_CFL_50_NX_201.png) |
| 3 | MacCormack ; supersonic outlet ; CFL = 1.1 ; Nx = 201 | not converged at final saved iteration | [P03](P03_MACCORMACK_SUPER_CFL_110_NX_201.png), [C03](C03_MACCORMACK_SUPER_CFL_110_NX_201.png) |
| 4 | MacCormack ; subsonic outlet ; CFL = 0.1 ; Nx = 201 | invalid flow at iteration 4498 | [P04](P04_MACCORMACK_SUB_CFL_10_NX_201.png), [C04](C04_MACCORMACK_SUB_CFL_10_NX_201.png) |
| 5 | MacCormack ; subsonic outlet ; CFL = 0.5 ; Nx = 201 | near stopping tolerance (CSV precision limited) | [P05](P05_MACCORMACK_SUB_CFL_50_NX_201.png), [C05](C05_MACCORMACK_SUB_CFL_50_NX_201.png) |
| 6 | MacCormack ; subsonic outlet ; CFL = 1.1 ; Nx = 201 | near stopping tolerance (CSV precision limited) | [P06](P06_MACCORMACK_SUB_CFL_110_NX_201.png), [C06](C06_MACCORMACK_SUB_CFL_110_NX_201.png) |
| 7 | Beam-Warming ; supersonic outlet ; CFL = 0.1 ; Nx = 201 | near stopping tolerance (CSV precision limited) | [P07](P07_BEAM_WARMING_SUPER_CFL_10_NX_201.png), [C07](C07_BEAM_WARMING_SUPER_CFL_10_NX_201.png) |
| 8 | Beam-Warming ; supersonic outlet ; CFL = 0.5 ; Nx = 201 | near stopping tolerance (CSV precision limited) | [P08](P08_BEAM_WARMING_SUPER_CFL_50_NX_201.png), [C08](C08_BEAM_WARMING_SUPER_CFL_50_NX_201.png) |
| 9 | Beam-Warming ; supersonic outlet ; CFL = 1.1 ; Nx = 201 | near stopping tolerance (CSV precision limited) | [P09](P09_BEAM_WARMING_SUPER_CFL_110_NX_201.png), [C09](C09_BEAM_WARMING_SUPER_CFL_110_NX_201.png) |
| 10 | Beam-Warming ; subsonic outlet ; CFL = 0.1 ; Nx = 201 | invalid flow at iteration 3893 | [P10](P10_BEAM_WARMING_SUB_CFL_10_NX_201.png), [C10](C10_BEAM_WARMING_SUB_CFL_10_NX_201.png) |
| 11 | Beam-Warming ; subsonic outlet ; CFL = 0.5 ; Nx = 201 | not converged at final saved iteration | [P11](P11_BEAM_WARMING_SUB_CFL_50_NX_201.png), [C11](C11_BEAM_WARMING_SUB_CFL_50_NX_201.png) |
| 12 | Beam-Warming ; subsonic outlet ; CFL = 1.1 ; Nx = 201 | not converged at final saved iteration | [P12](P12_BEAM_WARMING_SUB_CFL_110_NX_201.png), [C12](C12_BEAM_WARMING_SUB_CFL_110_NX_201.png) |
| 13 | MacCormack ; supersonic outlet ; CFL = 0.5 ; Nx = 101 | near stopping tolerance (CSV precision limited) | [P13](P13_MACCORMACK_SUPER_CFL_50_NX_101.png), [C13](C13_MACCORMACK_SUPER_CFL_50_NX_101.png) |
| 14 | MacCormack ; supersonic outlet ; CFL = 0.5 ; Nx = 401 | near stopping tolerance (CSV precision limited) | [P14](P14_MACCORMACK_SUPER_CFL_50_NX_401.png), [C14](C14_MACCORMACK_SUPER_CFL_50_NX_401.png) |
| 15 | MacCormack ; subsonic outlet ; CFL = 0.5 ; Nx = 101 | near stopping tolerance (CSV precision limited) | [P15](P15_MACCORMACK_SUB_CFL_50_NX_101.png), [C15](C15_MACCORMACK_SUB_CFL_50_NX_101.png) |
| 16 | MacCormack ; subsonic outlet ; CFL = 0.5 ; Nx = 401 | not converged at final saved iteration | [P16](P16_MACCORMACK_SUB_CFL_50_NX_401.png), [C16](C16_MACCORMACK_SUB_CFL_50_NX_401.png) |
| 17 | Beam-Warming ; supersonic outlet ; CFL = 0.5 ; Nx = 101 | near stopping tolerance (CSV precision limited) | [P17](P17_BEAM_WARMING_SUPER_CFL_50_NX_101.png), [C17](C17_BEAM_WARMING_SUPER_CFL_50_NX_101.png) |
| 18 | Beam-Warming ; supersonic outlet ; CFL = 0.5 ; Nx = 401 | near stopping tolerance (CSV precision limited) | [P18](P18_BEAM_WARMING_SUPER_CFL_50_NX_401.png), [C18](C18_BEAM_WARMING_SUPER_CFL_50_NX_401.png) |
| 19 | Beam-Warming ; subsonic outlet ; CFL = 0.5 ; Nx = 101 | near stopping tolerance (CSV precision limited) | [P19](P19_BEAM_WARMING_SUB_CFL_50_NX_101.png), [C19](C19_BEAM_WARMING_SUB_CFL_50_NX_101.png) |
| 20 | Beam-Warming ; subsonic outlet ; CFL = 0.5 ; Nx = 401 | not converged at final saved iteration | [P20](P20_BEAM_WARMING_SUB_CFL_50_NX_401.png), [C20](C20_BEAM_WARMING_SUB_CFL_50_NX_401.png) |

## Comparisons

- G01, CFL, runs [1, 2, 3]: [profiles](G01_CFL_profiles.png), [convergence](G01_CFL_convergence.png)
- G02, CFL, runs [4, 5, 6]: [profiles](G02_CFL_profiles.png), [convergence](G02_CFL_convergence.png)
- G03, CFL, runs [7, 8, 9]: [profiles](G03_CFL_profiles.png), [convergence](G03_CFL_convergence.png)
- G04, CFL, runs [10, 11, 12]: [profiles](G04_CFL_profiles.png), [convergence](G04_CFL_convergence.png)
- G05, grid, runs [13, 2, 14]: [profiles](G05_grid_profiles.png), [convergence](G05_grid_convergence.png)
- G06, grid, runs [15, 5, 16]: [profiles](G06_grid_profiles.png), [convergence](G06_grid_convergence.png)
- G07, grid, runs [17, 8, 18]: [profiles](G07_grid_profiles.png), [convergence](G07_grid_convergence.png)
- G08, grid, runs [19, 11, 20]: [profiles](G08_grid_profiles.png), [convergence](G08_grid_convergence.png)
