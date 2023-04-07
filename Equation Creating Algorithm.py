import sympy as sy
import math
import numpy as np

def Rot(φ):
  rad = math.radians(φ)
  return np.array([[math.cos(rad), -math.sin(rad)], [math.sin(rad), math.cos(rad)]])
#회전변환 행렬

A = float(input("x 반지름 :"))
B = float(input("y 반지름 :"))
mx = float(input("중심 x 좌표 :"))
my = float(input("중심 y 좌표 :"))
φ = float(input("회전각(°) :"))

x, y = sy.symbols("x y")
#x, y 문자로 연산하기로 선언

a = Rot(φ) @ np.array([[x - mx], [y-my]])
#선형 변환

E1 = sy.Eq(((a[0][0])**2)/(A**2) + ((a[1][0])**2)/(B**2), 1)
#타원방정식

sy.plot_implicit(E1, (x, -10, 10), (y, -10, 10))
#음함수 그래프 출력 (출력 범위 : (-10,10), (-10, 10))