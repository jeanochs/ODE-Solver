# Unit and Integration Tests

## Unit Tests

### Checking Solver Setup Functions

Main Functions/Subroutines Tested:

1. Element Creation Subroutines
    1. L2 Element: Nodes 1 and 2: x = 1, 3
        **Shape Functions and Derivatives**
        $ N_0(-0.3) = 0.65 $
        $ N_1(0.43) = 0.715 $
        $ N_0'(-0.1) = -0.5 $
        $ N_1'(-0.5) = 0.5 $

        **Jacobian**
        $ J = (3 - 1)/2 = 1 $
    
        **Isoparametric Relation**
        $ x(0.1) = 2.1 $

    2. L3 Element: Nodes 1, 2, and 3: x = 0, 1.3, 3
        **Shape Functions and Derivatives**
        $ N_0(-0.1) = 0.055 $
        $ N_1(-0.2) = 0.96 $
        $ N_2(0.33) = 0.2194 $
        $ N_0'(0.13) = -0.37 $
        $ N_1'(0.5) = -1 $
        $ N_2'(0.6) = 1.1 $

        **Jacobian**
        $ J(0.15) = 1.56 $
    
        **Isoparametric Relation**
        $ x(0.2) = 1.608 $
   
2. Calculation Routines 
    Driving Function Used: $ f(x) = x^2 + x + 3 $
    ODE Coefficients Used: $ a = 4, b = 4 $
    1. L2 Element: Nodes 1 and 2: x = 1, 3
        **Composition Functions**
        $ \Pi_{00}^K(-0.75) = 25/16 $
        $ \Pi_{10}^K(0.75) = -25/16 $
        $ \Pi_{01}^K(0.1) = 41/25 $
        $ \Pi_{11}^K(-0.2) = 163/100 $

        $ \Pi_{0}^F(-0.35) = 4.97644 $
        $ \Pi_{1}^F(0.2) = 6.024 $

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
        $ \Pi_{00}^K(-0.1) = 0.1322 $
        $ \Pi_{10}^K(0.3) = -1.273 $
        $ \Pi_{20}^K(0) = -0.1667 $
        $ \Pi_{01}^K(-0.5) = 2.1933 $
        $ \Pi_{11}^K(0.65) = 0.3051 $
        $ \Pi_{21}^K(0.9) = -6.3022 $
        $ \Pi_{02}^K(-0.75) = -0.6911 $
        $ \Pi_{12}^K(-0.1) = 1.3786 $
        $ \Pi_{22}^K(-0.53) = 0.0956 $

        $ \Pi_{0}^F(-0.3) = 1.2436 $
        $ \Pi_{1}^F(0.35) = 11.9016 $
        $ \Pi_{1}^F(0.55) = 7.3036 $

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


### Checking Parser Functions

Unit Tests: 5

Main Functions/Subroutines Tested:
`parse_input_file`

1. `linear_mesh_1,in`: Passing Test

Check:
    Node Count = 4
    Element Count = 3
    Calculations:
        For Element 1 (Index 0):
            Node IDs 1 and 2 (Indices 0 and 1)
            $ N_0(0) = \frac{1}{2} $ and $ N_1(0) = \frac{1}{2} $
            $ N_0'(0) = -\frac{1}{2} $ and $ N_1'(0) = \frac{1}{2} $
            $ J(0) = \frac{1 - 0}{2} = \frac{1}{2} $
        For Element 3 (index 2):
            Node IDs 3 and 4 (Indices 2 and 3)
            $ N_0(0) = \frac{1}{2} $ and $ N_1(0) = \frac{1}{2} $
            $ N_0'(0) = -\frac{1}{2} $ and $ N_1'(0) = \frac{1}{2} $
            $ J(0) = \frac{3 - 2}{2} = \frac{1}{2} $

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


    
