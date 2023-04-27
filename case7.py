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

def angle(a, x, y, φ):
    a1 = Rot(-φ) @ np.array([[a[0] - x], [a[1] - y]])
    if math.atan2(a1[1], a1[0]) < 0:
        return math.atan2(a1[1], a1[0]) + 2 * math.pi
    else:
        return math.atan2(a1[1], a1[0])
# 교점 a=(a1,a2)와 중심의 x좌표, y좌표 반시계방향으로 기울어진 정도를 넣어 기준 축에서 교점까지의 각을 구할 수 있다.

def f(θ, A, B):
    if θ == math.pi / 2 or θ == 3 * math.pi / 2:
        return θ
    else:
        return math.acos(A * math.tan(θ) / B)
    
def ESA(θ_1, θ_2, A, B):
    return A * B * (f(θ_2, A, B) - f(θ_1, A, B))/ 2

def tri(A, B, C):
    return abs(np.linalg.det(np.array[[A[0], A[1], 1], [B[1],B[2],1], [C[1],C[2],1]])) / 2
#세 점의 좌표로 삼각형의 넓이 구하기

def area_sec_1(a, b):
    return ESA(angle(Eq_real[a],mx_1,my_1,φ_1), angle(Eq_real[b],mx_1,my_1,φ_1), A_1, B_1)
def area_sec_2(a, b):
    return ESA(angle(Eq_real[a],mx_2,my_2,φ_2), angle(Eq_real[b],mx_2,my_2,φ_2), A_2, B_2)
#두 교점 사이의 부채꼴 넓이 구하기
    
def area_seg_1(a, b):
    return area_sec_1(a, b) - tri(Eq_real[a], Eq_real[b], (mx_1, my_1))
def area_seg_2(a, b):
    return area_sec_2(a, b) - tri(Eq_real[a], Eq_real[b], (mx_2, my_2))
#두 교점 사이의 활꼴 넓이 구하기
    
if len(Eq_real) == 3:
    Result = min(area_seg_1(0,1), area_seg_2(0,1)) + min(area_seg_1(1,2), area_seg_2(1,2)) + min(area_seg_1(2,0), area_seg_2(2,0)) + tri(Eq_real[0],Eq_real[1], Eq_real[2])
    print(Result)
else: print("error")
