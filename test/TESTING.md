# Unit and Integration Tests

## Preface

For this document, the composition functions are the integrands to the integrals for the coefficient matrix and constant vector.

$$ \Pi_{ij}^K = \int_{-1}^1 {\phi_i'(\zeta) \phi_j'(\zeta){\cdot}\frac{1}{J} + a\phi_i(\zeta) \phi_j'(\zeta) + b\phi_i(\zeta) \phi_j(\zeta){\cdot}J \ d\zeta} $$

$$ \Pi_i^F = \int_{-1}^1 {f(x(\zeta)){\cdot}\phi_i(\zeta){\cdot}J \ d\zeta} $$

## Unit Tests

### Checking Solver Setup Functions

Main Functions/Subroutines Tested:

1. Element Creation Subroutines

    1. L2 Element: Nodes 1 and 2: x = 1, 3

        **Shape Functions and Derivatives**

        $$ N_0(-0.3) = 0.65 $$ 

        $$ N_1(0.43) = 0.715 $$

        $$ N_0'(-0.1) = -0.5 $$

        $$ N_1'(-0.5) = 0.5 $$

        **Jacobian**

        $$ J = (3 - 1)/2 = 1 $$
    
        **Isoparametric Relation**

        $$ x(0.1) = 2.1 $$

    2. L3 Element: Nodes 1, 2, and 3: x = 0, 1.3, 3

        **Shape Functions and Derivatives**

        $$ N_0(-0.1) = 0.055 $$

        $$ N_1(-0.2) = 0.96 $$

        $$ N_2(0.33) = 0.2194 $$

        $$ N_0'(0.13) = -0.37 $$

        $$ N_1'(0.5) = -1 $$

        $$ N_2'(0.6) = 1.1 $$

        **Jacobian**
        $$ J(0.15) = 1.56 $$
    
        **Isoparametric Relation**
        $$ x(0.2) = 1.608 $$
   
2. Calculation Routines 

    Driving Function Used: $$ f(x) = x^2 + x + 3 $$

    ODE Coefficients Used: $$ a = 4, b = 4 $$

    1. L2 Element: Nodes 1 and 2: x = 1, 3

        **Composition Functions**

        $$ \Pi_{00}^K(-0.75) = 25/16 $$

        $$ \Pi_{10}^K(0.75) = -25/16 $$

        $$ \Pi_{01}^K(0.1) = 41/25 $$

        $$ \Pi_{11}^K(-0.2) = 163/100 $$

        $$ \Pi_{0}^F(-0.35) = 4.97644 $$

        $$ \Pi_{1}^F(0.2) = 6.024 $$

        **Vector and Matrix Outputs**
        $$ K = \begin{bmatrix}
            7/6 & 17/6 \\
            -7/6 & 31/6 \\
            \end{bmatrix} $$

        $$ F = \begin{bmatrix}
            23/3 \\
            11 \\
            \end{bmatrix} $$

    2. L3 Element: Nodes 1, 2 and 3: x = 0, 1.3, 3

        **Composition Functions**
        $$ \Pi_{00}^K(-0.1) = 0.1322 $$

        $$ \Pi_{10}^K(0.3) = -1.273 $$

        $$ \Pi_{20}^K(0) = -0.1667 $$

        $$ \Pi_{01}^K(-0.5) = 2.1933 $$

        $$ \Pi_{11}^K(0.65) = 0.3051 $$

        $$ \Pi_{21}^K(0.9) = -6.3022 $$

        $$ \Pi_{02}^K(-0.75) = -0.6911 $$

        $$ \Pi_{12}^K(-0.1) = 1.3786 $$

        $$ \Pi_{22}^K(-0.53) = 0.0956 $$

        $$ \Pi_{0}^F(-0.3) = 1.2436 $$

        $$ \Pi_{1}^F(0.35) = 11.9016 $$

        $$ \Pi_{2}^F(0.55) = 7.3036 $$

        **Vector and Matrix Outputs**
        $$ K = \begin{bmatrix}
            0.2099 & 2.2006 & -0.9438 \\
            -3.133 & 8.2577 & 2.8750 \\
            0.3895 & -2.4583 & 4.6021 \\
            \end{bmatrix} $$

        $$ F = \begin{bmatrix}
            0.4423 \\
            13.7783 \\
            8.2794 \\
            \end{bmatrix} $$

### Function Field Checks

Using $ f(x) = 3x^2 - 5x + 3 $ as the generating function, from 0 to 10, by steps 0.1, 101 number of points:

$$ f(2) = 5 $$

$$ f(5.5) = 66.25 $$

$$ f(4.67) \approx 45.0767 $$

## Integration Tests

### Parser Checks

#### Linear Meshes

1. `linear_mesh_1,in`: Passing Test

Check:
    Node Count = 4
    Element Count = 3
    Calculations:
        For Element 1 (Index 0):
            Node IDs 1 and 2 (coordinates x = 0 and x = 1)
            $$ N_0(0) = \frac{1}{2} $$ and $$ N_1(0) = \frac{1}{2} $$
            $$ N_0'(0) = -\frac{1}{2} $$ and $$ N_1'(0) = \frac{1}{2} $$
            $$ J(0) = \frac{1 - 0}{2} = \frac{1}{2} $$
        For Element 3 (index 2):
            Node IDs 3 and 4 (Coordinates x = 2 and x = 3)
            $$ N_0(0) = \frac{1}{2} $$ and $$ N_1(0) = \frac{1}{2} $$
            $$ N_0'(0) = -\frac{1}{2} $$ and $$ N_1'(0) = \frac{1}{2} $$
            $$ J(0) = \frac{3 - 2}{2} = \frac{1}{2} $$

2. `linear_mesh_2.in` : Failing Test

Check:
    Should error out; node count doesn't match.
    Return error code 1.

3. `linear_mesh_3.in` : Failing Test

Check:
    Should error out; nodes 3 and 4 (indices 2 and 3) are out of order.
    Return error code 1.

4. `linear_mesh_4.in` : Failing Test

Check:
    Should error out; should report that only the node count is given, and no node coordinates.
    Return error code 1.

5. `linear_mesh_5.in` : Failing Test

Check:
    Should error out; should report that file is empty.
    Return error code 1.

#### Quadratic Meshes

1. `quadratic_mesh_1.in` : Passing Test

Check:
    Node Count = 9
    Element Count = 4
    Calculations:
        For Element 1 (Index 0):
            Node IDs 1, 2 and 3 (coordinates x = 0, x = 1.25, and x = 2.5)
            $$ N_0(0.1) = -0.045 $$, $$ N_1(-0.3) = 0.91 $$, $$ N_2(0.5) = 0.375 $$
            $$ N_0'(0.73) = 0.23 $$, $$ N_1'(-0.34) = 0.68 $$, $$ N_2'(0.01) = 0.51 $$
            $$ J(-0.5) = 1.25  $$
            $$ x(0.4) = 1.75 $$
        For Element 3 (index 2):
            Node IDs 3 and 4 (Coordinates x = 2 and x = 3)
            $$ N_0(-0.43) = 0.30745 $$, $$ N_1(-0.99) = 0.02 $$, $$ N_2(0.48) = 0.3552 $$
            $$ N_0'(0.13) = -0.37 $$, $$ N_1'(-0.13) = 0.26 $$, $$ N_2'(0.57) = 1.07 $$
            $$ J(-0.07) = 0.542 $$
            $$ x(0.1) = 6.847 $$

2. `quadratic_mesh_2.in` : Failing Test

Check:
    Should error out; node count doesn't match.
    Return error code 1.

3. `quadratic_mesh_3.in` : Failing Test

Check:
    Should error out; nodes 3 and 4 (indices 2 and 3) are out of order.
    Return error code 1.

4. `quadratic_mesh_4.in` : Failing Test

Check:
    Should error out; should report that only the node count is given, and no node coordinates.
    Return error code 1.

5. `quadratic_mesh_5.in` : Failing Test

Check:
    Should error out; should report that file is empty.
    Return error code 1.

### Solver Checks

For both cases, the following mathematical parameters are used:

**Driving Function**

$$ f(x) = x^2 + x + 3 $$

**Constants**

$$ a = 4 $$

$$ b = 4 $$

$$ d_1 = 0 $$

$$ d_2 = 5 $$

1. Linear Mesh Solution

The coefficient matrix and constant vectors are as follows:

**Coefficient Matrix**

$$
\left[\begin{matrix}0.63188 & 3.9681 & 0 & 0 & 0\\-0.031884 & 2.1366 & 3.8952 & 0 & 0\\0 & -0.10476 & 2.4381 & 4.0667 & 0\\0 & 0 & 0.066667 & 8.7111 & 5.2222\\0 & 0 & 0 & 1.2222 & 7.7778\end{matrix}\right]
$$

**Constant Matrix**

$$
\left[\begin{matrix}5.3456\\12.447\\32.051\\172.52\\190.97\end{matrix}\right]
$$

The solution vector is given as:

**Solution Vector**

$$
\left[\begin{matrix}0\\30.869\\-13.737\\16.912\\5.0\end{matrix}\right]
$$

2. Quadratic Mesh Solution 

**Coefficient Matrix**

$$
\left[\begin{matrix}-1.6 & 4.4 & -1.1333 & 0 & 0 & 0 & 0 & 0 & 0\\-0.93333 & 3.2 & 4.4 & 0 & 0 & 0 & 0 & 0 & 0\\0.2 & -0.93333 & -5.2378 & 8.3552 & -1.7174 & 0 & 0 & 0 & 0\\0 & 0 & 3.0219 & 0.36339 & 5.9481 & 0 & 0 & 0 & 0\\0 & 0 & -0.38407 & 0.61475 & 4.2406 & -0.703 & 3.3651 & 0 & 0\\0 & 0 & 0 & 0 & -6.0363 & 21.915 & -13.212 & 0 & 0\\0 & 0 & 0 & 0 & 4.6984 & -18.545 & 12.536 & 4.3556 & -1.1778\\0 & 0 & 0 & 0 & 0 & 0 & -0.97778 & 4.6222 & 4.3556\\0 & 0 & 0 & 0 & 0 & 0 & 0.15556 & -0.97778 & 2.8222\end{matrix}\right]
$$

**Constant Matrix**

$$
\left[\begin{matrix}0.98958\\10.208\\-2.0435\\50.242\\66.22\\35.667\\27.6\\168.4\\56.05\end{matrix}\right]
$$

The solution vector is given as:

**Solution Vector**

$$
\left[\begin{matrix}0\\69.276\\-48.062\\-15.684\\33.822\\-5.112\\-26.632\\26.087\\5.0\end{matrix}\right]
$$





