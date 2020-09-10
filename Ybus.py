import numpy as np

z12 = complex(0.05, 0.25)
z23 = complex(0.02, 0.15)
z25 = complex(0.025, 0.1)
z34 = complex(0.01, 0.1)
z35 = complex(0.02, 0.1)

y12 = 1/z12
y23 = 1/z23
y25 = 1/z25
y34 = 1/z34
y35 = 1/z35

y11 = y12
y22 = y12 + y23 + y25
y33 = y34 + y23 + y35
y44 = y34
y55 = y25 + y35


y_bus = np.matrix([[y11, -y12, 0, 0, 0], [-y12, y22, -y23, 0, -y25], [0, -y23, y33, -y34, -y35], [0, 0, -y34, y44, 0],
                [0, -y25, -y35, 0, y55]])
