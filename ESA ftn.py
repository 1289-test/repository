#defining Area of Ellipse Sector for given θ_1, θ_2, major_semi_axis, minor_semi_axis

import math

def f(θ, A, B):
    if θ == math.pi / 2 or θ == 3 * math.pi / 2:
        return θ
    else:
        return math.acos(A * math.tan(θ) / B)
    
def ESA(θ_1, θ_2, A, B):
    return A * B * (f(θ_2, A, B) - f(θ_1, A, B))/ 2

    