import biomc_pp
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import svd
from scipy.optimize import minimize

np.set_printoptions(precision=1, suppress=True)
ds = "meta_batch_0d"


def fix_zero_rows(P):
    N = P.shape[0]
    for i in range(N):
        if P[i].sum() == 0:
            P[i] = np.ones(N) / N
    return P


# def infer_markov_matrix_from_histograms(h_t, h_t1, band_width=None):

#     n = len(h_t)

#     h_t=h_t/np.sum(h_t)
#     h_t1=h_t1/np.sum(h_t1)

#     P_init = np.random.rand(n, n)
#     P_init /= P_init.sum(axis=1, keepdims=True)

#     def loss(P_flat):
#         P = P_flat.reshape(n, n)
#         return np.sum((h_t1 - P @ h_t) ** 2)

#     def row_sum_constraint(P_flat):
#         P = P_flat.reshape(n, n)
#         return np.sum(P, axis=1) - 1

#     constraints = [{'type': 'eq', 'fun': row_sum_constraint}]
#     bounds = [(0, 1)] * (n * n)

#     res = minimize(loss, P_init.flatten(), bounds=bounds, constraints=constraints, method='SLSQP')

#     if not res.success:
#         print("failed")
#             # raise RuntimeError("Optimization failed: " + res.message)

#     P_t = res.x.reshape(n, n)
#     return P_t
#


def infer_markov_matrix_from_histograms(h_t, h_t1, band_width=2):
    n = len(h_t)
    h_t = h_t / np.sum(h_t)
    h_t1 = h_t1 / np.sum(h_t1)

    # P_init = np.eye(n)
    # P_init = np.zeros((n, n))
    # for i in range(n):
    #     for j in range(max(0, i - band_width), min(n, i + band_width + 1)):
    #         P_init[i, j] = np.random.rand()
    #
    P_init = np.random.random((n, n))

    P_init /= P_init.sum(axis=1, keepdims=True)

    def loss(P_flat):
        P = P_flat.reshape(n, n)
        return np.sum((h_t1 - P @ h_t) ** 2)

    def row_sum_constraint(P_flat):
        P = P_flat.reshape(n, n)
        return np.sum(P, axis=1) - 1

    constraints = [{"type": "eq", "fun": row_sum_constraint}]
    epsilon = 1e-6
    bounds = [(0, 1)] * (n * n)
    # bounds = []
    # for i in range(n):
    #     for j in range(n):
    #         if i+2<j:
    #             bounds.append((0, epsilon))
    #         else:
    #             bounds.append((0, 1))

    res = minimize(
        loss, P_init.flatten(), bounds=bounds, constraints=constraints, method="SLSQP"
    )

    if not res.success:
        raise RuntimeError("Optimization failed: " + res.message)

    P_t = res.x.reshape(n, n)

    return P_t


def _svd(K, bins):
    U, S, V_T = svd(K)
    plt.figure(figsize=(6, 4))
    plt.plot(S, "o-", label="Singular Value")
    plt.yscale("log")
    plt.xlabel("Mode index")
    plt.ylabel("Singular Valuee")
    plt.title("SVD Spectrum")
    plt.legend()
    plt.grid()
    lastmode = 3
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for i in range(3):
        uval = i * int((lastmode) / (3))

        axes[i].plot(bins[:-1], U[:, uval], "o-", label=f"Mode {i+1} (U)")
        axes[i].set_title(f"Mode {uval} - Initial impact")
        axes[i].set_xlabel("Length")
        axes[i].set_ylabel("Magnitude")
        axes[i].legend()
        axes[i].grid()

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for i in range(3):
        vval = i * int((lastmode) / (3))
        axes[i].plot(bins[:-1], V_T[vval, :], "o-", label=f"Mode {i+1} (Vt)")
        axes[i].set_title(f"Mode {vval} - Final Initial impact")
        axes[i].set_xlabel("Length")
        axes[i].set_ylabel("Magnitude")
        axes[i].legend()
        axes[i].grid()

    # # fig, ax = plt.subplots(4, 1, figsize=(8, 20))
    # S = np.diag(S)
    # for r in [5, 10, 70, 100]:
    #     fig = plt.figure()
    #     approx = U[:, :r] @ S[0:r, :r] @ V_T[:r, :]
    #     plt.imshow(approx, cmap="hot", aspect="auto")
    #     # ax[curr_fig].set_title("k = "+str(r))
    #     # ax[curr_fig].axis('off')

    # plt.show()


def dist_to_transition(histograms):
    nt, num_bins = histograms.shape
    T = np.zeros((num_bins, num_bins))

    histograms = histograms / np.sum(histograms, axis=1, keepdims=True)

    for t in range(nt - 1):
        for i in range(num_bins):
            if histograms[t, i] > 0:
                T[i, :] += histograms[t, i] * histograms[t + 1, :]

    T = T / np.sum(T, axis=1, keepdims=True)

    return [fix_zero_rows(np.nan_to_num(T, nan=0))]


pp = biomc_pp.PostProcess(ds, "/home-local/casale/Documents/code/poc/results/")
# def transition_matrix_sequence(H):
#     nt, N = H.shape
#     P_t_list = []

#     for t in range(nt - 1):
#         P = np.nan_to_num(H/H.sum(axis=1, keepdims=True),nan=0)

#         P_t_list.append(fix_zero_rows(P))


#         return P_t_list
#
def get_histogram(nt, nbins, lrange=None):
    H = np.zeros((nt, nbins))
    hist, bin_edges = np.histogram(
        pp.get_properties("length", 0) * 1e6, bins=nbins, range=lrange, density=False
    )
    H[0, :] = hist
    j = 0
    for i in range(0, nt):
        x = pp.get_properties("length", i) * 1e6
        hist, bin_edges = np.histogram(x, bins=bin_edges, range=lrange, density=False)
        H[j, :] = hist
        j += 1
    return H, bin_edges


def histogram_to_count_matrix(H):
    nt, nbins = H.shape
    count_matrices = []

    for t in range(nt - 1):
        h_t = H[t, :]
        h_next = H[t + 1, :]
        M = np.zeros((nbins, nbins))
        total_next = h_next.sum()
        if total_next > 0:
            for i in range(nbins):
                if h_t[i] > 0:
                    M[i, :] = h_t[i] * (h_next / total_next)
        count_matrices.append(M)
    return count_matrices


# def histogram_to_count_matrix(H):
#     nt, nbins = H.shape
#     count_matrices = []
#     for t in range(nt - 1):
#         h_t = H[t, :]
#         h_next = H[t+1, :]
#         M = np.zeros((nbins, nbins), dtype=int)
#         for i in range(nbins):
#             for j in range(nbins):
#                M[i, j] = min(h_t[i], h_next[j])
#         count_matrices.append(M)
#     return count_matrices


def count_to_transition_matrix(N):
    nbins = N.shape[0]
    P = np.zeros_like(N)
    for i in range(nbins):
        row_sum = N[i, :].sum()
        if row_sum != 0:
            P[i, :] = N[i, :] / row_sum
        else:
            P[i, :] = np.ones(nbins) / nbins
    return P


def transition_matrix_sequence_from_counts(count_matrices):
    P_list = []
    for N in count_matrices:
        P = count_to_transition_matrix(N)
        P_list.append(P)
    return P_list


def transition_matrix_sequence_from_hist(H):
    nt, nbins = H.shape
    P_list = []
    histograms = H / np.sum(H, axis=1, keepdims=True)
    for t in range(nt - 1):
        P = np.zeros((nbins, nbins))
        h_t = H[t, :]
        h_next = H[t + 1, :]
        row_sum = h_t.sum()

        if row_sum != 0:
            for i in range(nbins):
                if h_t[i] != 0:
                    P[i, :] = h_next / h_t[i]
                else:
                    P[i, :] = 1
            for i in range(nbins):
                row_sum_P = P[i, :].sum()
                if row_sum_P != 0:
                    P[i, :] /= row_sum_P
        else:
            P[:] = np.ones((nbins, nbins)) / nbins

        # T = np.zeros((nbins, nbins))
        # for i in range(nbins):
        #     if histograms[t, i] > 0:
        #         T[i, :] += histograms[t, i] * histograms[t + 1, :]
        # sum_T = np.sum(T, axis=1, keepdims=True)
        # mask_sum = np.where(sum_T!=0)
        # T[mask_sum]=T[mask_sum]/sum_T[mask_sum]
        # T = T / np.sum(T, axis=1, keepdims=True)

        P_list.append(P)
    return P_list


def markov_chain(P_list, n1, n2):
    K = P_list[n1]
    for i in range(n1, n2 - 1):
        K = K @ P_list[i]
    return K


lrange = (1.5, 3)
nbins = 50
nt = pp.n_export  # pp.n_export#np.maximum(pp.n_export,10)
H, bin_edges = get_histogram(nt, nbins, lrange)
Pt = dist_to_transition(H)  # transition_matrix_sequence_from_hist(H) #
# K = Pt[-1]#markov_chain(Pt,0,nt)#
# Pt = [
#   infer_markov_matrix_from_histograms(H[i, :], H[i + 1, :])
#   for i in range(450, nt - 1)
# ]
# # K = infer_markov_matrix_from_histograms(H[0,:],H[1,:])
# K = markov_chain(Pt, 0, len(Pt))  #
# print(np.sum(K, axis=1))

# plt.figure(figsize=(6, 5))
# plt.imshow(
#     K, cmap="hot", aspect="auto", extent=(lrange[0], lrange[1], lrange[1], lrange[0])
# )
# plt.colorbar(label="Proba")
# plt.ylabel("Init ")
# plt.xlabel("Final")
# plt.title("Markov Transition")


plt.phase_spectrum(np.cos(biomc_pp.check_time_unit(pp)))
plt.show()

#plt.figure(figsize=(6, 5))
#plt.imshow(
#     Pt[0],
#     cmap="hot",
#     aspect="auto",
#     extent=(lrange[0], lrange[1], lrange[1], lrange[0]),
# )
# plt.colorbar(label="Proba")
# plt.ylabel("Init ")
# plt.xlabel("Final")
# plt.title("Markov Transition 0 ")


# plt.figure(figsize=(6, 5))
# plt.imshow(Pt[-1], aspect="auto", extent=(lrange[0], lrange[1], lrange[1], lrange[0]))
# plt.colorbar(label="Proba")
# plt.ylabel("Init ")
# plt.xlabel("Final")
# plt.title("Markov Transition -1")


# l0 = H[0,:]
# estimation = l0
# jl = nt-1
# for i in range(1, nt, int((nt+1)/10)):
#     print(estimation,H[i,:])
#     estimation = estimation @ Pt[i]
# # estimation = l0@K
# #
# #
# fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# ax[0].bar(bin_edges[:-1], H[jl,:], alpha=0.3, label="H ")
# ax[0].bar(bin_edges[:-1], estimation, alpha=0.2, label="Estimation", color='r')
# ax[0].legend()
# ax[1].bar(bin_edges[:-1], (estimation - H[jl,:]) ** 2, label="Error", color='g')
# ax[1].legend()

# plt.tight_layout()
# plt.show()


_svd(K,bin_edges)
plt.figure()
for i in range(1, nt, int((nt + 1) / 10)):
    plt.hist(
        pp.get_properties("length", i) * 1e6, bins=100, density=False, label=f"{i}"
    )
plt.legend()
plt.figure()
# for i in range(1, nt):
#     l =pp.get_properties("length",i)*1e6
#     plt.hist(l,bins=100,density=False,label=f"{i}")


plt.show()
# plt.figure()
# pp = biomc_pp.PostProcess("udf_batch_2","/home-local/casale/Documents/code/poc/results/")
# # plt.hist(pp.get_properties("age",pp.n_export-1),bins=100,label="10",density=False)
# for i in range(1, pp.n_export, int((pp.n_export+1)/10)):
#     plt.hist(pp.get_properties("length",i)*1e6,bins=100,density=False,label=f"{i}")
# plt.title("batch")
# plt.legend()


# plt.figure()
# plt.plot(pp.get_number_particle())


# def gaussian(x, mean=0, std=1):
#     return np.exp(-0.5 * ((x - mean) / std) ** 2) / (std * np.sqrt(2 * np.pi))

# def dirac(x, epsilon=1e-3):
#     return np.exp(-0.5 * (x / epsilon) ** 2) / (epsilon * np.sqrt(2 * np.pi))

# def rect(x,midthd):
#     return np.where(abs(x)<=midthd, 1, 0)

# n = 1000
# x = np.linspace(0, 5, n)

# gaussian_signal = gaussian(x, mean=2.5, std=0.1)
# dirac_signal = dirac(x-3.5/2)

# result = convolve(gaussian_signal, dirac_signal, mode='same', method='auto')


# plt.figure(figsize=(10, 6))
# plt.subplot(3, 1, 1)
# plt.plot(x, gaussian_signal, label='Gaussian')
# plt.subplot(3, 1, 2)
# plt.plot(x, dirac_signal, label='Dirac')
# plt.subplot(3, 1, 3)
# plt.plot(x, result, label='Convolution')
# plt.tight_layout()
