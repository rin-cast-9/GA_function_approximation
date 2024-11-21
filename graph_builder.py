import matplotlib.pyplot as plt
import numpy as np

def read_data(file_path):
    with open(file_path, 'r') as file:
        step = float(file.readline().strip())

        y_target = list(map(float, file.readline().strip().split()))
        y_best = list(map(float, file.readline().strip().split()))

    return step, y_target, y_best


def plot_graphs(step, y_target, y_best):
    x_values = np.arange(0, 10 + step, step)

    plt.figure(figsize=(10, 5))
    plt.plot(x_values, y_target, label="target", linestyle='-')
    plt.plot(x_values, y_best, label="best solution", linestyle='-')

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)
    plt.show()


def main():
    # file_path = "graph_data.txt"
    file_path = "/Users/vadimabramov/Uni_works/programming/ai_systems/genetic_algorithm/graph_data.txt"
    plot_graphs(* read_data(file_path))


if __name__ == "__main__":
    main()
