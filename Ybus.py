import numpy as np

z12 = complex(0.05, 0.2)
z13 = complex(0.05, 0.1)
z23 = complex(0.05, 0.15)


y12 = 1/z12
y13 = 1/z13
y23 = 1/z23


y11 = y12 + y13
y22 = y12 + y23
y33 = y13 + y23


y_bus = np.matrix([[y11, -y12, -y13], [-y12, y22, -y23], [-y13, -y23, y33]])

#print(y_bus)