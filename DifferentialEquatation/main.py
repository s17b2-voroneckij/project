import matplotlib.pyplot as plt
import numpy as np
from bisect import bisect_left


def f_der2(x, t):
    return 2 * (x ** 2) + 2


def solve(n, a, b, x_0, f_der):  # f(a) = x_0, f_der = df/dt
    h = (b - a) / n
    x = [x_0] + [0] * n
    for k in range(1, n + 1):
        x[k] = x[k - 1] + h * f_der(x[k - 1], a + h * (k - 1))
    return [a + i * h for i in range(n + 1)], x


def pr7_a():
    eps = 0.1 / 2
    t1 = np.arange(0, np.pi / 8 - eps, 0.001)
    x1 = np.tan(2 * t1 + np.pi / 4)
    print(x1[0], np.pi / 8)
    colors = ['blue', 'red', 'black', 'yellow', 'purple', 'orange', 'green', 'pink', 'violet', 'brown']
    num_steps = [10, 30, 50, 100, 300, 500, 1000]
    for i in range(len(num_steps)):
        t2, x2 = solve(num_steps[i], 0, np.pi / 8 - eps, x1[0], f_der2)
        plt.plot(t2, x2, color=colors[i + 1])

    plt.plot(t1, x1, color=colors[0])
    plt.show()


def draw_slopes(f_der, a, C, de):
    for dt in range(-a, a + 1):
        for dx in np.arange(-a, a + 1, 1):
            t = dt / de
            x = dx / de
            der = f_der(x, t)
            if der != np.inf:
                T = np.array([-C, C])
                X = np.array([-C * der, C * der])
                len = np.sqrt((T[0] - T[1]) ** 2 + (X[0] - X[1]) ** 2) / C
                T /= len
                X /= len
                T += t
                X += x
                plt.plot(T, X, color='red')
            else:
                plt.plot([t, t], [x - C / 2, x + C / 2], color='red')


def draw_vector_field(f_x, f_y, a, d):
    for dx in range(-a, a + 1):
        for dy in range(-a, a + 1):
            x = dx / d
            y = dy / d
            plt.plot([x, x + f_x(x)], [y, y + f_y(y)], color='red')


def make_plot(f, a, b):
    t = np.arange(a, b, 0.01)
    x = f(t)
    plt.plot(t, x, color='blue')


def f1(x, t): return 1 / 3 * t / x if x != 0 else np.inf


def f2(x, t): return -5 * t / x if x != 0 else np.inf


def f3(x, t): return -1 / 4 * x / t if t != 0 else np.inf


def f4(x, t): return 1 / 5 * x / t if t != 0 else np.inf


def f6_1(x, t): return (9 * t - 3 * x - 5) ** 2 + 2


def f6_2(x, t): return 3 - 3 * (x ** 2)


def pr1():
    draw_slopes(f1, 30, 0.05, 8)
    make_plot(lambda t: np.sqrt((t ** 2) / 3 + 4), -4, 4)
    make_plot(lambda t: -np.sqrt((t ** 2) / 3 + 1), -4, 4)
    plt.show()


def pr2():
    draw_slopes(f2, 30, 0.05, 8)
    make_plot(lambda t: np.sqrt(4 - 5 * (t ** 2)), -2 / np.sqrt(5), 2 / np.sqrt(5))
    make_plot(lambda t: -np.sqrt(1 - 5 * (t ** 2)), -1 / np.sqrt(5), 1 / np.sqrt(5))
    plt.show()


def pr3():
    draw_slopes(f3, 30, 0.05, 8)
    make_plot(lambda t: 1 / abs(t) ** (1 / 4), -4, -0.01)
    make_plot(lambda t: 1 / abs(t) ** (1 / 4), 0.01, 4)
    make_plot(lambda t: -1 / abs(t) ** (1 / 4), -4, -0.01)
    make_plot(lambda t: -1 / abs(t) ** (1 / 4), 0.01, 4)
    make_plot(lambda t: t - t, -4, 4)
    plt.show()


def pr4():
    draw_slopes(f4, 30, 0.05, 8)
    make_plot(lambda t: abs(t) ** (1/5), -4, 0)
    make_plot(lambda t: -abs(t) ** (1 / 5), -4, 0)
    make_plot(lambda t: abs(t) ** (1 / 5), 0, 4)
    make_plot(lambda t: -abs(t) ** (1 / 5), 0, 4)
    plt.show()


def pr5_a():
    draw_vector_field(lambda x: -x, lambda y: 4 * y, 5, 2)
    make_plot(lambda x: 1 / (x ** 4), -4, -0.6)
    make_plot(lambda x: 1 / (x ** 4), 0.6, 4)
    plt.plot([0] * 2, [-12, 12], color='black')
    plt.show()


def pr5_b():
    draw_vector_field(lambda x: x, lambda y: 5 * y, 5, 2)
    make_plot(lambda x: x ** 5, -2, 2)
    plt.plot([0] * 2, [-12, 12], color='black')
    plt.show()


def pr6_2():
    draw_slopes(f6_2, 30, 0.05, 8)

    plt.show()


def pr6_1():
    draw_slopes(f6_1, 30, 0.05, 8)
    # make_plot(lambda t: 3 * t - (np.e ** (6 * t) - 1) / 3 / (np.e ** (6 * t) + 1) - 5 / 3, -1, 2)
    # make_plot(lambda t: 3 * t - 2, -1, 2)
    # make_plot(lambda t: 3 * t - 4 / 3, -1, 2)
    make_plot(lambda t: 3 * t - (40 + np.e ** (6 * t)) / 3 / (np.e ** (6 * t) - 40) - 5 / 3, -1, -0.1)
    for st in [-300, -250, -200, -150, -50, -40, -10, 1, 2, 5, 10, 50, 100, 200, 250]:
        make_plot(lambda t: 3 * t - (10 ** (st) + np.e ** (6 * t)) / 3 / (np.e ** (6 * t) - 10 ** (st)) - 5 / 3, -1, -0.1)
    plt.show()


def get_value(point, T, X):
    i = bisect_left(T, point)
    return X[i]


if __name__ == '__main__':
    pr6_1()
