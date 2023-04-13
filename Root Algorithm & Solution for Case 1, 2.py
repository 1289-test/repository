import sympy as sy
import math
import numpy as np

def Rot(φ):
  rad = math.radians(φ)
  return np.array([[math.cos(rad), -math.sin(rad)], [math.sin(rad), math.cos(rad)]])
#회전변환 행렬

A_1 = float(input("x 반지름 :"))
B_1 = float(input("y 반지름 :"))
mx_1 = float(input("중심 x 좌표 :"))
my_1 = float(input("중심 y 좌표 :"))
φ_1 = float(input("(시계 방향)회전각(°) :"))

x, y = sy.symbols("x y")
#x, y 를 연산 가능한 기호로 선언

a = Rot(φ_1) @ np.array([[x - mx_1], [y-my_1]])
#변환

E1 = sy.Eq(((a[0][0])**2)/(A_1**2) + ((a[1][0])**2)/(B_1**2) - 1 , 0)
#타원방정식 a는 행렬꼴, matrix[n+1번째 entry][0 = eq type] 

A_2 = float(input("x 반지름 :"))
B_2 = float(input("y 반지름 :"))
mx_2 = float(input("중심 x 좌표 :"))
my_2 = float(input("중심 y 좌표 :"))
φ_2 = float(input("(시계 방향)회전각(°) :"))

b = Rot(φ_2) @ np.array([[x - mx_2], [y-my_2]])
#변환

E2 = sy.Eq(((b[0][0])**2)/(A_2**2) + ((b[1][0])**2)/(B_2**2) - 1, 0)
#타원방정식

Eq = sy.solve((E1, E2), (x, y))
#두 타원의 교점 (x,y) 꼴 딕셔너리

Eq_real = []
#일반적으로 앞서 x, y의 수 체계를 default 로 설정함
#sympy 모듈에서의 기호의 default 수 체계로 복소수를 취하므로 복소근을 삭제해야함 
#실근 딕셔너리를 따로 만듦

for i in range(len(Eq)):
    if Eq[i][0].is_real:
        Eq_real.append(Eq[i])
    else:
        pass
#전체 근 딕셔너리에서 x값이 실수이면 실근 딕셔너리에 추가

#Eq_real, len(Eq_real)
#이로써 실근의 좌표와 개수를 알 수 있음.

if len(Eq_real) = 0:
    #####CASE 1########################################################################
    if E1.sy.subs([x, mx_2], [y, my_2]) < 0:#E2의 중점을 E1에 대입했을때, 내부에 있는가?
        Result = A_2 * B_2 * math.pi() / 2 #E1이 E2를 포함하는 형태 즉 E2의 넓이를 return
    elif E2.sy.subs([x, mx_1], [y, mx_2]) < 0:#반대
        Result = A_1 * B_1 * math.pi / 2
    #####CASE 2########################################################################
    else:#포함관계가 없는 경우
        Result = 0
##### 이후론 CASE 3~#####################################################################
elif len(Eq_real) = 1:
    Result = 'none' #CASE 3, 4...
#... 