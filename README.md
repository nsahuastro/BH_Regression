# Black Hole Mass Regression Analysis

This project analyzes the relationship between galaxy properties and supermassive black hole masses using statistical regression techniques. 

The code Mbh_M_ETG_LTG.py focuses on two main galaxy types: Early-Type Galaxies (ETGs) and Late-Type Galaxies (LTGs).

## What does the script do?

- Loads galaxy data from a CSV file.
- Filters and separates the data into ETGs and LTGs.
- Performs regression analysis to find correlations between galaxy mass and black hole mass.
- Calculates statistical metrics such as scatter, intrinsic scatter, and correlation coefficients.
- Plots the results, showing regression lines and confidence intervals for both galaxy types.

## Who is this for?

This script is intended for astronomers, astrophysicists, and students interested in galaxy evolution and black hole scaling relations. Basic familiarity with Python and scientific data analysis is helpful.

## How to use

1. Place your galaxy data CSV file (e.g., `New_ETG_LTG3.csv`) in the project directory.
2. Ensure required Python packages are installed: `numpy`, `pandas`, `matplotlib`, `scipy`, and any custom modules referenced (e.g., `bces`, `stats`, `misc`, `xplot`).
3. Run the script `Mbh_M_ETG_LTG.py` using Python 3:
   ```
   python Mbh_M_ETG_LTG.py
   ```
4. The script will print statistical results and display a plot comparing ETGs and LTGs.

## Output

- Statistical summaries for regression fits.
- A plot showing the relationship between galaxy mass and black hole mass for ETGs and LTGs.

## Notes

- Some galaxies are excluded from analysis based on specific criteria.
- The script uses advanced regression methods suitable for astronomical data.
- Modify the CSV filename or filtering criteria as needed for your dataset.

## Contact

For questions or suggestions, please contact the project maintainer.

## References

- Sahu, N., Graham, A. W., Davis, B. L., et al. (2019a). [The Supermassive Black Hole Mass–Galaxy Mass Relation for Early- and Late-Type Galaxies](https://iopscience.iop.org/article/10.3847/1538-4357/ab0f32/meta#apjab0f32s4). *The Astrophysical Journal*, 876, 155.
- [BCES Regression Python Implementation](https://github.com/rsnemmen/BCES)


