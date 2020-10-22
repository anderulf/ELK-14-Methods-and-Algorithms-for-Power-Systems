import numpy as np
#Assignment 4 part 1
B_p=np.array ([[15, -5, -10],[ -5, 9, -4],[-10,-4, 18]])
X_12 = 0.2
X_13 = 0.1
X_23 = 0.25
X_34 = 0.25

P_1 = -1.25
P_2 = -0.4
P_3 = -0.6

#bus 4 is slack!

# PTDF


# By assuming no change in voltages, one implies zero reactive power flow

# line 1-2
b=np.array([[1/X_12],[-1/X_12],[0]])
a=np.matmul(np.linalg.inv(B_p),b)
P=a[0,0]*P_1+a[1,0]*P_2+a[2,0]*P_3

P_change1 = a[0,0]*(P_1-0.5)+a[1,0]*P_2+a[2,0]*P_3
P_change2 = a[0,0]*(P_1-0.5)+a[1,0]*(P_2+0.3)+a[2,0]*P_3
print('-----Line 1-2 ----')
print("b:",b)
print("a:",a)
print("P_12:",P)
print('When increasing the load at bus one with 0.5pu')
print("The loadflow 1-2 is canged by:",P_change1-P)
print('When increasing the load at bus one with 0.5pu, and decreased by 0.3 at bus 2')
print("The loadflow 1-2 is canged by:",P_change2-P)
print(' ')
# line 1-3
b=np.array([[1/X_13],[0],[-1/X_13]])
a=np.matmul(np.linalg.inv(B_p),b)
P=a[0,0]*P_1+a[1,0]*P_2+a[2,0]*P_3
P_change1 = a[0,0]*(P_1-0.5)+a[1,0]*P_2+a[2,0]*P_3
P_change2 = a[0,0]*(P_1-0.5)+a[1,0]*(P_2+0.3)+a[2,0]*P_3
print('-----Line 1-3 ------')
print("b:",b)
print("a:",a)
print("P_13:",P)
print('When increasing the load at bus one with 0.5pu')
print("The loadflow 1-3 is canged by:",P_change1-P)
print('When increasing the load at bus one with 0.5pu, and decreased by 0.3 at bus 2')
print("The loadflow 1-3 is canged by:",P_change2-P)

print(' ')
# line 2-3
b = np.array([[0],[1/X_23],[-1/X_23]])
a = np.matmul(np.linalg.inv(B_p),b)
P = a[0,0]*P_1+a[1,0]*P_2+a[2,0]*P_3

P_change1 = a[0,0]*(P_1-0.5)+a[1,0]*P_2+a[2,0]*P_3
P_change2 = a[0,0]*(P_1-0.5)+a[1,0]*(P_2+0.3)+a[2,0]*P_3
print('------Line 2-3 -------')
print("b:",b)
print("a:",a)
print("P_23:",P)
print('When increasing the load at bus one with 0.5pu')
print("The loadflow 2-3 is canged by:",P_change1-P)
print('When increasing the load at bus one with 0.5pu, and decreased by 0.3 at bus 2')
print("The loadflow 2-3 is canged by:",P_change2-P)
print(' ')

# line 3-4
b = np.array([[0],[0],[-1/X_34]])
a = np.matmul(np.linalg.inv(B_p),b)
P = a[0,0]*P_1+a[1,0]*P_2+a[2,0]*P_3

P_change1 = a[0,0]*(P_1-0.5)+a[1,0]*P_2+a[2,0]*P_3
P_change2 = a[0,0]*(P_1-0.5)+a[1,0]*(P_2+0.3)+a[2,0]*P_3
print('-----Line 3-4------ ')
print("b:",b)
print("a:",a)
print("P_34:",P)


print('\nThe consept of distribution factors is to avoid calculating the angles')
print('look above to see the distribution factor for each line (a)')
print('When increasing the load at bus one with 0.5pu')
print("The loadflow on line 3-4 is changed by: {}pu".format(round(P_change1-P, 3)) )
print('When increasing the load at bus one with 0.5pu, and decreasing by 0.3pu at bus 2')
print("The loadflow on line 3-4 is changed by: {}pu".format(round(P_change2-P, 3)))

print(' ')

print("\n Task 1")
