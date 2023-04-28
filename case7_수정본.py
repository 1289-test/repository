import sympy as sy
import math
import numpy as np
import decimal as dec

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
E_1 = [A_1, B_1, mx_1, my_1, φ_1]
E1_deriv = sy.idiff(E1.lhs, y, x)

#타원방정식 a는 행렬꼴, matrix[n+1번째 entry][0 = eq type] 

A_2 = float(input("x 반지름 :"))
B_2 = float(input("y 반지름 :"))
mx_2 = float(input("중심 x 좌표 :"))
my_2 = float(input("중심 y 좌표 :"))
φ_2 = float(input("(시계 방향)회전각(°) :"))

b = Rot(φ_2) @ np.array([[x - mx_2], [y-my_2]])
#변환

E2 = sy.Eq(((b[0][0])**2)/(A_2**2) + ((b[1][0])**2)/(B_2**2) - 1, 0)
E_2 = [A_2, B_2, mx_2, my_2, φ_2]
E2_deriv = sy.idiff(E2.lhs, y, x)
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

def angle(a, x, y, φ):
    a1 = Rot(-φ) @ np.array([[a[0] - x], [a[1] - y]])
    if math.atan2(a1[1], a1[0]) < 0:
        return math.atan2(a1[1], a1[0]) + 2 * math.pi
    else:
        return math.atan2(a1[1], a1[0])
# 교점 a=(a1,a2)와 중심의 x좌표, y좌표 반시계방향으로 기울어진 정도를 넣어 기준 축에서 교점까지의 각을 구할 수 있다.

def sort_angle(Eq_list):
    for n in range(len(Eq_list)-1):
        for m in range(n + 1, len(Eq_list)):
            if angle(Eq_list[n], mx_1, my_1, φ_1) > angle(Eq_list[m], mx_1, my_1, φ_1):
                Eq_list[n], Eq_list[m] = Eq_list[m], Eq_list[n]           
#Eq_real을 각도 순으로 정렬 (첫번째 타원 기준)

def f(θ, A, B):
    if θ == math.pi / 2 or θ == 3 * math.pi / 2:
        return θ
    elif 0 <= θ < math.pi / 2:
        return math.atan(A * math.tan(θ) / B)
    elif math.pi / 2 < θ <= math.pi:
        return math.atan(A * math.tan(θ) / B) + math.pi
    elif math.pi < θ < 3 * math.pi / 2:
        return math.atan(A * math.tan(θ) / B) + math.pi
    elif 3 * math.pi / 2 < θ < 2 * math.pi:
        return math.atan(A * math.tan(θ) / B) + 2 * math.pi
    
def ESA(θ_1, θ_2, A, B):
    return A * B * (f(θ_2, A, B) - f(θ_1, A, B))/ 2

def area_tri(a, b, mx, my):
    return abs((a[0] - mx) * (b[1] - my) - (b[0] - mx) * (a[1] - my)) / 2    
# 두 교점과 타원의 중심의 좌표를 넣으면 세 점이 이루는 삼각형의 넓이를 반환함

def area_seg(a, b, E):
    if ESA(angle(a, E[2], E[3], E[4]), angle(b, E[2], E[3], E[4]), E[0], E[1]) < 0:
        if E[0] * E[1] * math.pi + ESA(angle(a, E[2], E[3], E[4]), angle(b, E[2], E[3], E[4]), E[0], E[1]) > E[0] * E[1] * math.pi / 2:
            return E[0] * E[1] * math.pi + ESA(angle(a, E[2], E[3], E[4]), angle(b, E[2], E[3], E[4]), E[0], E[1]) + area_tri(a, b, E[2], E[3])
        elif E[0] * E[1] * math.pi + ESA(angle(a, E[2], E[3], E[4]), angle(b, E[2], E[3], E[4]), E[0], E[1]) < E[0] * E[1] * math.pi / 2:
            return E[0] * E[1] * math.pi + ESA(angle(a, E[2], E[3], E[4]), angle(b, E[2], E[3], E[4]), E[0], E[1]) - area_tri(a, b, E[2], E[3])
        else: 
            return E[0] * E[1] * math.pi / 2
    elif ESA(angle(a, E[2], E[3], E[4]), angle(b, E[2], E[3], E[4]), E[0], E[1]) > E[0] * E[1] * math.pi / 2:
        return ESA(angle(a, E[2], E[3], E[4]), angle(b, E[2], E[3], E[4]), E[0], E[1]) + area_tri(a, b, E[2], E[3])
    elif ESA(angle(a, E[2], E[3], E[4]), angle(b, E[2], E[3], E[4]), E[0], E[1]) < E[0] * E[1] * math.pi / 2:
        return ESA(angle(a, E[2], E[3], E[4]), angle(b, E[2], E[3], E[4]), E[0], E[1]) - area_tri(a, b, E[2], E[3])
    else:
        return E[0] * E[1] * math.pi / 2
#두 교점의 좌표와 타원리스트를 넣으면 두 교점이 만드는 활꼴의 넓이를 반환함. (첫번째 교점에서 두번째 교점으로 반시계방향으로)
    
if len(Eq_real) == 3:
    ####CASE 7###
    E=sort_angle(Eq_real)
    Result = min(area_seg(E[0],E[1],E_1), area_seg(E[0],E[1],E_2)) +  min(area_seg(E[1],E[2],E_1), area_seg(E[1],E[2],E_2))\
        + min(area_seg(E[2],E[0],E_1), area_seg(E[2],E[0],E_2)) +  abs((E[0][0] - E[2][0]) * (E[1][1] - E[2][1]) - (E[1][0] - E[2][0]) * (E[0][1] - E[2][1])) / 2

    
    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    