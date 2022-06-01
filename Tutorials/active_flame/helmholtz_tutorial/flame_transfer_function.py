import numpy as np

def n_tau(N3, tau):
    """
    :param N3: non-dimensional interaction index
    :param tau: time delay [s]
    :return: function
    """
    def inner_func(omega, k=0):
        return N3 * (1j * tau)**k * np.exp(1j * omega * tau)
    return inner_func


# if __name__ == '__main__':
#     f=n_tau(0.16, 0.50)
#     print(f(0))
