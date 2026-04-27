from sage.all import *
import time
from tqdm import tqdm
import itertools
import pickle
from coxeter3_sage import Coxeter3


class Alg:
    # this class collect useful methods provide by SageMath and Coxeter3
    # At the moment only integral weights/level is implemented
    def __init__(self, cartanType, QLoad=True, WLoad=True):
        print(">>> Initialing Alg object <<<", flush=True)
        self.cartanType = cartanType
        # self.r denotes the finite algebra rank
        if len(cartanType) == 3:
            # self.rank is the affine rank
            self.rank = self.cartanType[1] + 1
            self.r = self.rank - 1
            self.weight_lattice = RootSystem(cartanType).weight_lattice(extended=True)
        else:
            self.rank = self.cartanType[1]
            self.r = self.rank
            self.weight_lattice = RootSystem(cartanType).weight_lattice()

        # basic quantities
        self.omega = self.weight_lattice.fundamental_weights()
        self.alpha = self.weight_lattice.simple_roots()
        self.rho = sum(self.omega)

        if len(cartanType) == 3:
            # finite highest root (θ; 0; 0)
            # translate back to an affine root in the second step
            self.theta = self.weight_lattice.classical().positive_roots_by_height()[-1]
            self.theta = sum(
                [self.omega[i] * self.theta.to_vector()[i - 1] for i in range(1, self.r + 1)]
            )

            # self.weight_lattice.basic_imaginary_roots()[0] is
            # the same as self.weight_lattice.null_root()
            self.delta = self.weight_lattice.basic_imaginary_roots()[0]
            # returns the classical Cartan matrix
            self.A = matrix(self.weight_lattice.root_system.cartan_matrix()[1:, 1:])
            if self.cartanType[0] == "G":
                self.A = self.A.transpose()
            self.finite_positive_roots = list(self.weight_lattice.classical().positive_roots())
            self.finite_fundamental_weights = self.weight_lattice.classical().fundamental_weights()
        elif len(cartanType) == 2:
            self.theta = self.weight_lattice.positive_roots_by_height()[-1]
            self.delta = 0
            self.A = matrix(self.weight_lattice.root_system.cartan_matrix())
            self.cartan_matrix = self.A
            self.finite_positive_roots = list(self.weight_lattice.positive_roots())
            self.finite_fundamental_weights = self.weight_lattice.fundamental_weights()
        self.AInversed = self.A.inverse()

        self.comarks = {}
        self.comarks.update(
            zip(
                [i for i in range(1, self.r + 1)],
                [self.ScalarProduct(self.theta, self.omega[i]) for i in range(1, self.r + 1)],
            )
        )
        self.hcheck = 1 + sum(self.comarks.values())

        # Weyl group implementations
        # self.WW uses the built-in weyl_group
        # self.W will be generated manually
        self.WW = self.weight_lattice.weyl_group(prefix="w")
        self.ws = self.WW.generators()  # simple reflections
        self.one = self.ws[0] * self.ws[0]  # identity

        # For affine Weyl group
        # prepare simple translations t_i = s_i * s_{αi + δ}
        if len(cartanType) == 3:
            print("\t>>> Getting simple translations <<<", flush=True)
            # w_{α[i] + δ}, i = 0, 1, ..., r = rank - 1
            # note that w_{α[0] + δ} is also present
            wiDelta = [
                self.WW.from_reduced_word((self.alpha[i] + self.delta).associated_reflection())
                for i in range(self.rank)
            ]
            self.simpleTranslations = [self.ws[i] * wiDelta[i] for i in range(self.rank)]

        # initialize self.W for manual generation through self.GetWeylGroup(depth)
        self.W = [self.one]
        self.T = [self.one]
        # self.WW.classical() are different objects from
        # those in self.WW
        # therefore can't perform product between the two types of elements
        # need to translate back to elements in self.WW exploiting the
        # reduced_word() expression
        print("\t>>> Getting finite Weyl group <<<", flush=True)
        start = time.time()
        if len(cartanType) == 3 and WLoad:
            try:
                print("\t\t>>> Loading finite Weyl group from file.", flush=True)
                self.load_W_finite()
            except:
                print("\t\t>>> File not found. Generating finite Weyl group.", flush=True)
                self.WFinite = [
                    self.WW.from_reduced_word(w.reduced_word()) for w in self.WW.classical().list()
                ]
                if len(self.WFinite) < 10000:
                    self.WFinite = sorted(
                        [
                            self.WW.from_reduced_word(w.reduced_word())
                            for w in self.WW.classical().list()
                        ],
                        key=lambda x: x.length(),
                    )
        elif len(cartanType) == 3 and not WLoad:
            print("\t\t>>> WLoad = False. Generating finite Weyl group.", flush=True)
            self.WFinite = [
                self.WW.from_reduced_word(w.reduced_word()) for w in self.WW.classical().list()
            ]
        elif len(cartanType) == 2:
            self.WFinite = self.WW
        end = time.time()
        print("\t>>> Took " + str(end - start) + "s", flush=True)
        print("\t>>> self.WFinite is set", flush=True)

        # some data for KL computation
        q = var("q")
        self.cox = Coxeter3(self.WW, q)
        try:
            self.WWCoxeter3 = CoxeterGroup(self.cartanType, implementation="coxeter3")
        except:
            pass
        self.llambda = 0 * self.rho
        self.Lambda = 0 * self.rho
        self.WLambda0 = []  # WLambda0 is used in Qtilde
        self.QData = {}
        # When an instance is created,
        # update the Qdata from reading the stored data
        # might take some time
        if QLoad == True:
            self.QLoad()
        print(">>> Initialization complete <<<", flush=True)

    def Dynkin(self, weight):
        # Returns the Dynkin labels of the finite part
        try:
            # for affine extended weight
            labels = weight.to_classical().to_vector()
            return list(labels)
        except:
            # for finite weight
            labels = weight.to_vector()
            return list(labels)

    def AffineDynkin(self, weight):
        # λ0 ωhat[0] + λi ωhat[i] + ... --> [λ0, λ1, ... ]
        # weight.to_vector() returns [λ0, λ1, ..., λr, n]
        # e.g., for affine su(2)
        # α0.to_vector() = (2,-2,1), α1.to_vector() = (-2, 2, 0)
        # α0^v.to_vector() = (1,0), α1^v.to_vector() = (0,1)
        # NOTE: [0:self.rank] ~ entry 0, 1, 2, ..., self.rank - 1
        return list(weight.to_vector()[0 : self.rank])

    def ScalarProduct(self, weight1, weight2):
        # compute scalar product of FINITE weights
        # if weight1, weight2 are affine weights (λ, k, n), (μ, k', n')
        # return ONLY (λ, μ)
        lambda1 = self.Dynkin(weight1)
        lambda2 = self.Dynkin(weight2)
        # r = rank of the finite part
        if len(self.cartanType) == 2:
            r = self.rank
        elif len(self.cartanType) == 3:
            r = self.rank - 1
        temp = sum(
            flatten(
                [
                    [
                        lambda1[i - 1]
                        * lambda2[j - 1]
                        * self.AInversed[i - 1, j - 1]
                        * (self.alpha[j].norm_squared())
                        / 2
                        for i in range(1, r + 1)
                    ]
                    for j in range(1, r + 1)
                ]
            )
        )
        if self.cartanType[0] == "G":
            temp = temp / 3
        return temp

    def ToAffineWeight(self, weight):
        # map ([λ1, λ2, ..., λr];k;n) to an affine weight
        # Only works for integral weights
        if len(weight) == 3 and len(weight[0]) == self.r:
            k = weight[1]
            n = int(weight[2])
            finite_part = sum([weight[0][i - 1] * self.omega[i] for i in range(1, self.r + 1)])
            lambda0 = int(k - self.ScalarProduct(finite_part, self.theta))
            return self.omega[0] * lambda0 + finite_part + n * self.delta
        elif len(weight) == 3 and len(weight[0]) != self.r:
            print("finite piece is not rank ", self.r, flush=True)
            return

    def Tolambdakn(self, weight):
        # λ0 ωhat[0] + λi ωhat[i] --> ([λ1, λ2, ...]; k; n)
        dynkinLabelsWithDelta = weight.to_vector()
        return (
            list(self.Dynkin(weight)),
            dynkinLabelsWithDelta[0] + self.ScalarProduct(weight, self.theta),
            dynkinLabelsWithDelta[-1],
        )

    def AffineScalarProduct(self, weight1, weight2):
        v1 = self.Tolambdakn(weight1)
        v2 = self.Tolambdakn(weight2)
        return self.ScalarProduct(weight1, weight2) + v1[2] * v2[1] + v1[1] * v2[2]

    def GetMoreWeylGroupElements(self):
        # Take the Cartesian product of the current self.W with all simple reflections
        # to get more Weyl group elements
        newW = set(list(map(product, list(itertools.product(self.ws, list(self.W))))))
        if newW != set(self.W):
            return self.W.union(newW)
        else:
            return self.W

    def GetMoreTranslations(self):
        newT = set(
            [
                product(tuple)
                for tuple in list(itertools.product(list(self.T), self.simpleTranslations))
            ]
        )
        self.T = set(self.T)
        if newT != self.T:
            self.T = self.T.union(newT)
            return self.T
        else:
            return self.T

    def GetTranslationsByLength(self, length):
        self.T = {self.one}
        for i in range(length):
            self.GetMoreTranslations()
        self.T = sorted(list(self.T), key=lambda x: x.length())
        return self.T

    def _CeilSqrtQQ(self, value):
        target = QQ(value)
        if target <= 0:
            return 0

        lower = 0
        upper = 1
        while QQ(upper) * QQ(upper) < target:
            lower = upper
            upper *= 2

        while lower + 1 < upper:
            mid = (lower + upper) // 2
            if QQ(mid) * QQ(mid) >= target:
                upper = mid
            else:
                lower = mid

        return upper

    def _FiniteCorootGramMatrix(self):
        indices = list(range(1, self.rank))
        gram = []
        for i in indices:
            row = []
            for j in indices:
                numerator = self.ScalarProduct(2 * self.alpha[i], 2 * self.alpha[j])
                denominator = self.ScalarProduct(self.alpha[i], self.alpha[i]) * self.ScalarProduct(
                    self.alpha[j], self.alpha[j]
                )
                row.append(QQ(numerator) / QQ(denominator))
            gram.append(row)
        return gram

    def _TranslationCoefficientRadius(self, level, linear_coeffs, gram, max_neg_shift):
        k = QQ(level)
        b = QQ(self._CeilSqrtQQ(sum(QQ(d) * QQ(d) for d in linear_coeffs)))
        G = matrix(QQ, gram)
        G_inv = G.inverse()
        frob_sq = sum(QQ(v) * QQ(v) for v in G_inv.list())
        lam_lower = QQ(1) / QQ(self._CeilSqrtQQ(frob_sq))

        if lam_lower <= 0:
            raise ValueError(
                "Failed to obtain a positive coercive bound for translation enumeration"
            )

        abs_k = abs(k)
        a = abs_k * lam_lower / QQ(2)
        if a <= 0:
            raise ValueError(
                "Failed to derive a positive quadratic bound for translation enumeration"
            )

        if k > 0:
            disc = b * b + QQ(4) * a * QQ(max_neg_shift)
            r_real = (b + QQ(self._CeilSqrtQQ(disc))) / (QQ(2) * a)
            return max(0, int(r_real) + 1)

        r_real = b / a
        return max(0, int(r_real) + 1)

    def _MinusDeltaN(self, level, linear_coeffs, gram, coeffs):
        m = [QQ(c) for c in coeffs]
        linear = sum(QQ(d) * mi for d, mi in zip(linear_coeffs, m))
        quadratic = QQ(0)
        for i, mi in enumerate(m):
            for j, mj in enumerate(m):
                quadratic += mi * QQ(gram[i][j]) * mj
        return linear + QQ(level) * quadratic / QQ(2)

    def _TranslationFromCoefficients(self, coeffs):
        factors = [
            self.simpleTranslations[index] ** int(coeff)
            for index, coeff in enumerate(coeffs, start=1)
            if int(coeff) != 0
        ]
        if not factors:
            return self.one
        return product(factors)

    def _TranslationsByNShiftBnB(self, weight, order, order_min=0, return_stats=False):
        max_neg_shift_value = QQ(order)
        min_neg_shift_value = QQ(order_min)
        if max_neg_shift_value < min_neg_shift_value:
            return {"translations": [], "stats": {}} if return_stats else []

        level = QQ(self.AffineDynkin(weight)[0] + self.ScalarProduct(weight, self.theta))
        if level == 0:
            raise ValueError("Translation enumeration by n-shift requires non-zero level")

        linear_coeffs = [QQ(value) for value in self.AffineDynkin(weight)[1 : self.rank]]
        gram = self._FiniteCorootGramMatrix()
        radius = self._TranslationCoefficientRadius(
            level=level,
            linear_coeffs=linear_coeffs,
            gram=gram,
            max_neg_shift=max_neg_shift_value,
        )

        n = len(linear_coeffs)
        if n == 0:
            translations = [self.one] if min_neg_shift_value <= 0 <= max_neg_shift_value else []
            if return_stats:
                return {"translations": translations, "stats": {"radius": 0}}
            return translations

        try:
            G = matrix(QQ, gram)
            d_vec = vector(QQ, linear_coeffs)
            center_vec = -(G.solve_right(d_vec)) / level
            center = [QQ(center_vec[i]) for i in range(n)]
        except Exception:
            center = [QQ(0) for _ in range(n)]

        coeffs = [0 for _ in range(n)]
        selected = {}
        stats = {
            "radius": int(radius),
            "dimension": int(n),
            "box_points": int((2 * int(radius) + 1) ** int(n)),
            "visited_leaves": 0,
            "pruned_branches": 0,
            "accepted": 0,
        }

        gram_qq = [[QQ(gram[i][j]) for j in range(n)] for i in range(n)]
        tail_quad_lower = [QQ(0) for _ in range(n + 1)]
        tail_quad_upper = [QQ(0) for _ in range(n + 1)]
        R = QQ(radius)
        for depth in range(n - 1, -1, -1):
            lower = tail_quad_lower[depth + 1]
            upper = tail_quad_upper[depth + 1]
            qii = level * gram_qq[depth][depth] / QQ(2)
            term = qii * R * R
            if term >= 0:
                upper += term
            else:
                lower += term
            for j in range(depth + 1, n):
                band = abs(level * gram_qq[depth][j]) * R * R
                lower -= band
                upper += band
            tail_quad_lower[depth] = lower
            tail_quad_upper[depth] = upper

        current_const = QQ(0)
        current_b = [QQ(v) for v in linear_coeffs]

        def partial_bounds(depth):
            if depth == n:
                return current_const, current_const

            lower = QQ(current_const) + tail_quad_lower[depth]
            upper = QQ(current_const) + tail_quad_upper[depth]
            for i in range(depth, n):
                delta = abs(current_b[i]) * R
                lower -= delta
                upper += delta
            return lower, upper

        ordered_values_by_dim = []
        base_values = list(range(-radius, radius + 1))
        for i in range(n):
            values = list(base_values)
            values.sort(key=lambda x: (abs(QQ(x) - center[i]), abs(x), x))
            ordered_values_by_dim.append(values)

        def dfs(depth, prefix_norm_sq):
            nonlocal current_const, current_b
            if prefix_norm_sq > radius * radius:
                stats["pruned_branches"] += 1
                return

            low, high = partial_bounds(depth)
            if high < min_neg_shift_value or low > max_neg_shift_value:
                stats["pruned_branches"] += 1
                return

            if depth == n:
                stats["visited_leaves"] += 1
                neg_shift = self._MinusDeltaN(level, linear_coeffs, gram, tuple(coeffs))
                if min_neg_shift_value <= neg_shift <= max_neg_shift_value:
                    key = tuple(int(c) for c in coeffs)
                    selected[key] = self._TranslationFromCoefficients(key)
                    stats["accepted"] += 1
                return

            for value in ordered_values_by_dim[depth]:
                coeffs[depth] = int(value)
                value_qq = QQ(value)
                old_const = current_const
                old_b = list(current_b)

                current_const = (
                    old_const
                    + current_b[depth] * value_qq
                    + level * gram_qq[depth][depth] * value_qq * value_qq / QQ(2)
                )
                for j in range(depth + 1, n):
                    current_b[j] = old_b[j] + level * gram_qq[j][depth] * value_qq

                dfs(depth + 1, prefix_norm_sq + int(value) * int(value))
                current_const = old_const
                current_b = old_b

            coeffs[depth] = 0

        dfs(0, 0)
        translations = sorted(selected.values(), key=lambda x: x.length())
        if return_stats:
            return {"translations": translations, "stats": stats}
        return translations

    def GetTranslationsBynShift(self, weight, order=1, max_m=5, order_min=None):
        print(">>> Getting translations for {} to order {}".format(weight, order), flush=True)
        self.order = order
        if order_min is None:
            order_min = 0
        max_neg_shift_value = QQ(order)
        min_neg_shift_value = QQ(order_min)
        if max_neg_shift_value < min_neg_shift_value:
            print(">>> Found translations =  []", flush=True)
            return []

        level = QQ(self.AffineDynkin(weight)[0] + self.ScalarProduct(weight, self.theta))
        linear_coeffs = [QQ(value) for value in self.AffineDynkin(weight)[1 : self.rank]]
        gram = self._FiniteCorootGramMatrix()
        radius = self._TranslationCoefficientRadius(
            level=level,
            linear_coeffs=linear_coeffs,
            gram=gram,
            max_neg_shift=max_neg_shift_value,
        )

        dimension = len(linear_coeffs)
        box_points = (2 * int(radius) + 1) ** int(dimension)
        if box_points > 2000000:
            result = self._TranslationsByNShiftBnB(
                weight,
                order=max_neg_shift_value,
                order_min=min_neg_shift_value,
                return_stats=False,
            )
            print(">>> Found translations = ", result, flush=True)
            return result

        ranges = [range(-radius, radius + 1) for _ in linear_coeffs]
        selected = {}
        for coeffs in itertools.product(*ranges):
            neg_shift = self._MinusDeltaN(level, linear_coeffs, gram, coeffs)
            if neg_shift < min_neg_shift_value or neg_shift > max_neg_shift_value:
                continue
            selected[tuple(int(coeff) for coeff in coeffs)] = self._TranslationFromCoefficients(
                coeffs
            )

        Ts = [selected[key] for key in sorted(selected.keys())]
        Ts = sorted(Ts, key=lambda x: x.length())
        print(">>> Found translations = ", Ts, flush=True)
        return Ts

    get_translations_by_n_shift = GetTranslationsBynShift

    def nShift(self, weight, m):
        # m = (0, m1, m2, ..., mr)
        # t1^m1 t2^m2 ... on (λ;k;n) shifts
        # n -> n + Δn
        # Δn = - Sum[m[i]λ[i],i] - (k/2) Sum[m[i]m[j](αv[i], αv[j]), i,j]
        # Be careful with the negative sign: when weight is dominant, Δn <= 0

        # self.ScalarProduct() just perform the finite ScalarProduct.
        # For an affine dominant weight, w1, ..., wr will always decrease n.
        alpha = self.alpha
        # k = λ[0] + (λ, θ)
        k = self.AffineDynkin(weight)[0] + self.ScalarProduct(weight, self.theta)
        print("k = {}".format(k))
        return -sum([m[i] * self.AffineDynkin(weight)[i] for i in range(1, self.rank)]) - (
            k / 2
        ) * sum(
            [
                sum(
                    [
                        m[i]
                        * m[j]
                        * self.ScalarProduct(2 * alpha[i], 2 * alpha[j])
                        * self.ScalarProduct(alpha[i], alpha[i]) ** (-1)
                        * self.ScalarProduct(alpha[j], alpha[j]) ** (-1)
                        for i in range(1, self.rank)
                    ]
                )
                for j in range(1, self.rank)
            ]
        )

    def GetWeylGroup(self, l, style="height"):
        # for finite Weyl group
        if len(self.cartanType) == 2:
            self.W = list(self.WW.list())
            return self.W

        # for affine Weyl group
        if style == "semi-direct-product":
            self.GetTranslationsByLength(l)
            self.W = [product(tuple) for tuple in itertools.product(self.WFinite, self.T)]
            self.W = sorted(self.W, key=lambda x: x.length())
            return self.W
        elif style == "qSeries":
            # manually set alg.T according to q level in character computation
            # independent of the l argument
            # the number of elements will be |W| x |self.T|
            self.W = [product(tuple) for tuple in itertools.product(self.WFinite, self.T)]
            if len(self.W) < 1000:
                self.W = sorted(self.W, key=lambda x: x.length())
            return self.W
        elif style == "height":
            self.W = flatten([list(self.WW.elements_of_length(i)) for i in range(l)])
            return self.W

    def GetWeylGroupForqSeries(self, weight=None, order=2, T=None):

        # Leave <weight> absent if alg.T is set manually, in this case
        # <order> is also not useful
        print(">>> Creating Weyl group for q-series expansion to order %s" % order, flush=True)
        start = time.time()
        if T is None:
            T = self.GetTranslationsBynShift(weight, order)
        # use existing alg.T to compute alg.W
        print("\t>>> Taking semi-direct product.", flush=True)
        start1 = time.time()
        W = [product(tuple) for tuple in itertools.product(self.WFinite, T)]
        end1 = time.time()
        print("\t>>> Semi-direct product completed: %s s" % str(end1 - start1), flush=True)
        end = time.time()
        print(
            ">>> Weyl group created: total %s s" % str(end - start),
            "\n>>> Group size: %s\n" % len(W),
            flush=True,
        )
        return W

    get_Weyl_group_for_q_series = GetWeylGroupForqSeries

    def SaveWFinite(self):
        with open("finite weyl group/%s.dat" % self.cartanType[:2], "wb") as file:
            pickle.dump([w.reduced_word() for w in self.WFinite], file)

    save_W_finite = SaveWFinite

    def LoadWFinite(self):
        file_name = "finite weyl group/%s.dat" % self.cartanType[:2]
        with open(file_name, "rb") as file:
            WFiniteReducedWord = pickle.load(file)
            self.WFinite = [self.WW.from_reduced_word(word) for word in WFiniteReducedWord]
            print("\t>>> File name = {}\n".format(file_name), flush=True)
            file.close()

    load_W_finite = LoadWFinite

    def SaveT(self, weight=None, order=None):
        if not weight:
            weight = self.llambda
        if not order:
            order = self.order
        print(">>> Saving translations for {} to order {}".format(weight, order), flush=True)
        HWString = str(self.llambda).replace("Lambda", "omega").replace("*", "")

        filename = "translations/{}-({})-qorder{}.dat".format(
            str(self.cartanType).replace(" ", ""), HWString, order
        )
        with open(filename, "wb") as file:
            pickle.dump([w.reduced_word() for w in self.T], file)
            print(">>> File name = {}".format(filename), flush=True)

    save_translations = SaveT

    def SaveW(self):
        HWString = str(self.llambda).replace("Lambda", "omega").replace("*", "")

        filename = "affine weyl group/{}-({})-qorder{}.dat".format(
            str(self.cartanType).replace(" ", ""), HWString, self.order
        )
        with open(filename, "wb") as file:
            pickle.dump([w.reduced_word() for w in self.W], file)

    save_W = SaveW

    def SaveWDen(self, order=None):
        if not order:
            order = self.order
        filename = "affine weyl group/{}-rho-qorder{}.dat".format(
            str(self.cartanType).replace(" ", ""), order
        )
        with open(filename, "wb") as file:
            pickle.dump([w.reduced_word() for w in self.W_denonminator], file)

    save_W_denominator = SaveWDen

    def SaveFinitePositiveRoots(self, order=None):
        print(">>> Saving finite positive roots to file.")
        # after calling self.Kazhdan_Lusztig_denominator, save the numerator info into a txt file
        roots_str = str(self.finite_positive_roots)
        roots_str = "{" + roots_str[1:]
        roots_str = roots_str[:-1] + "}"
        filename = "finite positive roots/" + str(self.cartanType).replace(" ", "") + ".txt"
        file2 = open(filename, "w")
        file2.write(roots_str)
        file2.close()
        print(">>> file name = ", filename, "\n")
        return None

    save_finite_positive_roots = SaveFinitePositiveRoots

    def LoadW(self, llambda, order):
        print(">>> Loading affine W from file.", flush=True)
        HWString = str(llambda).replace("Lambda", "omega").replace("*", "")

        filename = "affine weyl group/{}-({})-qorder{}.dat".format(
            str(self.cartanType).replace(" ", ""), HWString, order
        )
        with open(filename, "rb") as file:
            WReducedWord = pickle.load(file)
            self.W = [self.WW.from_reduced_word(word) for word in WReducedWord]
            file.close()
        print(">>> Successfully loaded affine Weyl from file.\n", flush=True)
        return self.W

    load_W = LoadW

    def LoadT(self, llambda, order):
        print(
            ">>> Loading translations for {} to order {} from file.".format(llambda, order),
            flush=True,
        )
        HWString = str(llambda).replace("Lambda", "omega").replace("*", "")

        filename = "translations/{}-({})-qorder{}.dat".format(
            str(self.cartanType).replace(" ", ""), HWString, order
        )
        with open(filename, "rb") as file:
            reduced_words = pickle.load(file)
            self.T = [self.WW.from_reduced_word(word) for word in reduced_words]
            file.close()
        print(">>> Successfully loaded translations from file.\n", flush=True)
        return self.T

    load_translations = LoadT

    def LoadWDen(self, order):
        print("\t>>> Loading W_denonminator from file.", flush=True)
        filename = "affine weyl group/{}-rho-qorder{}.dat".format(
            str(self.cartanType).replace(" ", ""), order
        )
        with open(filename, "rb") as file:
            WReducedWord = pickle.load(file)
            self.W_denonminator = [self.WW.from_reduced_word(word) for word in WReducedWord]
            file.close()
        print("\t>>> Successfully loaded affine Weyl_denominator from file.", flush=True)

    load_W_denominator = LoadWDen

    def RemoveDelta(self, weight):
        return sum(
            [
                list(weight.to_vector())[0 : self.rank][i] * list(self.omega)[i]
                for i in range(self.rank)
            ]
        )

    def Reflection(self, root):
        # returns the reflection operation associated to finite/affine root
        return self.WW.from_reduced_word(root.associated_reflection())

    def WeylToList(self, w):
        return w.reduced_word()

    def ExtractLastw(self, w):
        ind = w.reduced_word()
        if len(ind) == 0:
            return self.one
        elif len(ind) == 1:
            return [self.one, w]

        return [self.WW.from_reduced_word(ind[:-1]), self.WW.from_reduced_word([ind[-1]])]

    # def MaxRep(self, w, subgroup):
    #     temp = w
    #     coset = [w * s for s in subgroup]
    #     for ws in coset:
    #         if ws.length() > w.length():
    #             temp = ws
    #     return temp
    # def MinRep(self, w, subgroup):
    #     temp = w
    #     coset = [w * s for s in subgroup]
    #     for ws in coset:
    #         if ws.length() < w.length():
    #             temp = ws
    #     return temp

    def MaxRep(self, w, subgroup):
        return max([w * s for s in subgroup], key=lambda w: w.length())

    def MinRep(self, w, subgroup):
        return min([w * s for s in subgroup], key=lambda w: w.length())

    def BruhatDescendents(self, w):
        # list all weyl group element that is bruhat-larger-than-or-equal-to w
        return [wp for wp in self.W if w.bruhat_le(wp) and not wp.bruhat_le(w)]

    def P(self, x, y):
        # invoke the KL polynomial from Coxeter3
        WW = self.WWCoxeter3
        if x == self.one and y == self.one:
            return WW.kazhdan_lusztig_polynomial([], [])
        elif x == self.one:
            return WW.kazhdan_lusztig_polynomial([], self.WeylToList(y))
        elif y == self.one:
            return WW.kazhdan_lusztig_polynomial(self.WeylToList(x), [])
        else:
            return WW.kazhdan_lusztig_polynomial(self.WeylToList(x), self.WeylToList(y))

    def Pslash(self, x, y):
        # extract the (l(y) - l(x) - 1)/2 order term in the P(x, y)
        poly = self.P(x, y)
        if poly in ZZ:
            if 0 == (y.length() - x.length() - 1) / 2:
                return poly
            else:
                return 0
        q = poly.variables()[0]
        if poly.degree(q) == (y.length() - x.length() - 1) / 2:
            # note that Coxeter3's P(x, y) is Dense univariate polynomials over Z
            # its .leading_coefficient() doesn't need an argument
            return (poly.leading_coefficient()) * q ** (poly.degree(q))
        else:
            return 0

    def Qslash(self, x, y):
        # Qslash(x, y) = Pslash(x, y), according to Vos and Driel
        return self.Pslash(x, y)

    @cached_method
    def Q(self, x, y, method=None):
        # Q is the inverse polynomial of the KL polynomial P

        # method = "coxeter3": force the program to get invpol(x,y) using coxeter3
        #       = "recursive": directly compute Q(x, y) using recursion relations
        #       = None: load saved values when possible, and get invpol(x,y) if necessary
        q = var("q")
        # Q(x, y) = 0 if not x <= y
        if not x.bruhat_le(y):
            return 0
        # Q(e, y) = 1, eq (A.16)
        if x == self.one:
            return 1
        # Q(x, x) = 1, below eq (A.11)
        if x == y:
            return 1
        # when stlye = "coxeter3", force the program to
        # compute the invpol using coxeter3
        # without reading the saved values
        if method == "coxeter3":
            print("Using cox.invpol(x,y) directly.", flush=True)
            return self.cox.invpol(x, y)

        if method == "recursive":
            print("Computing Q by recursion; (x, y) = ", (x, y), flush=True)
            [yy, s] = self.ExtractLastw(y)
            if self.lt(x * s, x):
                # c = 1
                result = (
                    self.Q(x * s, yy)
                    - q * self.Q(x, yy)
                    + q
                    * sum(
                        [
                            self.Qslash(x, z) * self.Q(z, yy)
                            for z in self.WW.bruhat_interval(x, y)  # x < z <= y
                            if self.lt(x, z) and self.lt(z, z * s)
                        ]  # zs > z
                    )
                )
            else:
                # c = 0
                result = self.Q(x, yy)
            return result

        # Default: when method = None
        # read saved data if available
        if (x, y) in self.QData.keys():
            return self.QData[(x, y)]
        # If no available saved Data, compute using coxeter3
        result = self.cox.invpol(x, y)
        if floor(time.time()) % 10 == 0:
            print("Random status check: computing (x, y) = ", (x, y), flush=True)

        # Simplify the result when possible
        if "full_simplify" in dir(result):
            result = result.full_simplify()

        # Save new result of Q(x,y) to self.QData
        if (x, y) not in self.QData.keys() and x.bruhat_le(y):
            self.QData[(x, y)] = result
        return result

    def QSave(self):
        print("Saving QData to file.", str(self.cartanType) + ".dat", flush=True)
        with open(str(self.cartanType) + ".dat", "wb") as f:
            pickle.dump(self.QData, f)
            print("QData saved: file name = ", f.name, "\n", flush=True)
        return

    def QLoad(self):
        start = time.time()
        print("\t>>> Loading inverse KL polynomials Q from file <<<", flush=True)
        try:
            with open(str(self.cartanType) + ".dat", "rb") as f:
                self.QData = pickle.load(f)
        except:
            end = time.time()
            print("\t>>> File not found <<<", flush=True)
            return
        end = time.time()
        print("\t>>> inverse KL polynomials loaded: %s s" % str(end - start), flush=True)
        return

    # can't cache this method, since subgroup is a list
    # and is not hashable
    def Qtilde(self, x, y, subgroup):
        xbar = self.MaxRep(x, subgroup)
        coset = [y * s for s in subgroup]
        result = sum([self.Q(xbar, z) * int((-1) ** (xbar.length() - z.length())) for z in coset])
        if "full_simplify" in dir(result):
            result = result.full_simplify()
        return result

    def lt(self, x, y):
        # check if x < y
        return x.bruhat_le(y) and not y.bruhat_le(x)

    def GetLambda(self, llambda):
        # lambda is a reserved word, use llambda instead
        # =====================================================================
        # When k + hcheck > 0
        # find the unique weight Lambda of the form "Λ = w(λ + ρ) - ρ"
        # such that Λ + ρ is dominant
        # i.e., Λ + ρ has NON-NEGATIVE (can be zero) Dynkin labels
        # =====================================================================
        start = time.time()
        print(">>> Creating Lambda.", flush=True)
        self.llambda = llambda
        finite_coefficients = list((llambda).to_vector()[0 : self.rank])
        if all(coeff > 0 for coeff in finite_coefficients):
            self.wToLambda = self.one
            self.wTollambda = self.one
            print("wTollambda = ", self.wTollambda, ", Λ = ", self.llambda, flush=True)
            return llambda

        print("Finding Λ for generic λ", flush=True)
        rho = self.rho
        lambda_plus_rho = llambda + rho
        checked = 0
        max_length = 0
        max_steps = 1000

        while checked < max_steps:
            elements_of_length = list(self.WW.elements_of_length(max_length))
            elements_of_length.sort(key=lambda w: tuple(int(i) for i in w.reduced_word()))
            for w_to_Lambda in elements_of_length:
                checked += 1
                acted_weight = w_to_Lambda.action(lambda_plus_rho) - rho
                reduced_coefficients = list(acted_weight.to_vector()[0 : self.rank])
                if all(coeff >= -1 for coeff in reduced_coefficients):
                    self.Lambda = acted_weight
                    self.wToLambda = w_to_Lambda
                    self.wTollambda = self.wToLambda.inverse()
                    print("wTollambda = ", self.wTollambda, ", Λ = ", self.Lambda, flush=True)
                    end = time.time()
                    print(">>> Lambda created: %s s" % str(end - start), "\n", flush=True)
                    return acted_weight
                if checked >= max_steps:
                    break
            max_length += 1

        raise Exception("Failed to find dominant Lambda by bounded affine Weyl search")

    # alias
    get_lambda = GetLambda

    def GetWLambda0(self, Lambda):
        # returns the isometric subgroup that fixes Λ
        # depends on the size of self.W
        # WLambda0 is used in Qtilde computation
        # and is crucial to get it right
        start = time.time()
        print(">>> Creating WLambda0.")
        self.WLambda0 = [w for w in self.W if w.action(Lambda + self.rho) - self.rho == Lambda]
        if not self.one in self.WLambda0:
            self.WLambda0 = self.WLambda0 + [self.one]
        end = time.time()
        print("W_Lambda0 = ", self.WLambda0)
        print(">>> WLambda0 created: %s s" % str(end - start), "\n")
        return self.WLambda0

    # alias
    get_W_Lambda_0 = GetWLambda0

    def CharacterNum(self, llambda, order=None, Lambda=None, wTollambda=None):
        # computes the numerator of the KL formula
        # which is a sum over affine Weyl (sub)group
        # or sum over a lots of weights
        print(">>> Computing Kazhdan-Lusztig numerator.")
        rho = self.rho

        self.llambda = llambda
        if Lambda is not None:
            self.Lambda = Lambda
            self.wTollambda = wTollambda
        else:
            self.Lambda = self.GetLambda(llambda)
        Lambda = self.Lambda
        self.order = order
        # If an <order> param is specified, regenerate the self.T
        # and self.W based on the shift in n-value of Λ+ρ up to
        # the <order> param
        try:
            self.T = self.load_translations(llambda, order)
        except:
            print("\t>>> File not found. Building translations from scratch.")
            order_min = (
                self.Tolambdakn(self.Lambda + self.rho)[-1] - self.Tolambdakn(self.llambda)[-1]
            )
            self.T = self.get_translations_by_n_shift(
                Lambda + rho, order_min + order, order_min=None
            )
            print("\t>>> self.T")
            print(self.T)
        self.W = self.GetWeylGroupForqSeries(order=order, T=self.T)
        # if order is not specified, assume self.W is manually set
        W = self.W

        # After Getting a suitable/large enough Weyl group
        # Computes the invariant subgroup W0_Λ
        self.GetWLambda0(Lambda)

        print("======== Long computation begins ========")

        # dotted orbit of Λ (not λ)
        # restricted to wTollambda <= w' where wTollambda(Λ+ρ)-ρ = λ
        # LambdaOrbitUnderWeylDot = [(w.action(self.Lambda + rho) - rho) for w # in self.W if self.wTollambda.bruhat_le(w)]
        start = time.time()
        print(">>> Creating LambdaOrbitUnderWeylDot")
        LambdaOrbitUnderWeylDot = [(w.action(self.Lambda + rho) - rho) for w in self.W]
        end = time.time()
        print(">>> LambdaOrbitUnderWeylDot created: %s s." % str(end - start), "\n")

        # Main objective: find the cosets
        # Strategy: look at the dot-action image set W.(Λ + rho)
        # w's having the same w.(Λ + ρ) are in the same cosets
        # for any coset, use the shorted representative
        # scan the image set W.(Λ + ρ)
        # - record unique element w.(Λ+ρ)
        # - record the corresponding w
        # - if there is a shorter w, replace the existing w
        weightsToBeSummed = []
        cosets = []
        start = time.time()
        print(">>> Creating cosets.")
        for i in tqdm(range(len(LambdaOrbitUnderWeylDot)), leave=True):
            weight = LambdaOrbitUnderWeylDot[i]
            w = W[i]
            if weight not in weightsToBeSummed:
                weightsToBeSummed.append(weight)
                cosets.append(w)
            else:
                ind = weightsToBeSummed.index(weight)
                if w.length() < cosets[ind].length():
                    # keep the shorter w
                    cosets[ind] = w
        end = time.time()
        print(">>> cosets created: %s s" % str(end - start))
        print("weightsToBeSummed length = ", len(weightsToBeSummed), "\n")

        print(">>> Creating WeylToBeSummed.")
        start = time.time()
        # ordering [w] <= [w'] is determined using the min rep in each coset
        # wTollambda is defined to be min rep
        # cosets is also constructed using min rep
        WeylToBeSummed = [wp for wp in cosets if self.wTollambda.bruhat_le(wp)]
        # store the results for debugging
        self.cosets = cosets
        self.weightsToBeSummed = weightsToBeSummed
        self.WeylToBeSummed = WeylToBeSummed
        end = time.time()
        print(">>> WeylToBeSummed created: %s s" % str(end - start))
        print("WeylToBeSummed length = ", len(WeylToBeSummed), "\n")
        print(">>> Creating numerator.")
        start = time.time()
        chunk_size = 20  # the size of each chunk

        def process_chunk(chunk):
            return [
                {
                    wp.action(self.Lambda + rho) - rho: self.Qtilde(
                        self.wTollambda, wp, self.WLambda0
                    )
                }
                for wp in chunk
            ]

        self.num = []  # list to collect results
        for i in tqdm(range(0, len(WeylToBeSummed), chunk_size), leave=True):
            chunk = WeylToBeSummed[i : i + chunk_size]  # get chunk WeylToBeSummed
            self.num = self.num + process_chunk(chunk)

        end = time.time()
        print(">>> numerator created: %s s\n" % str(end - start))
        return self.num

    Kazhdan_Lusztig_numerator = CharacterNum

    def CharacterDen(self, order):
        print(">>> Creating denominator.", flush=True)
        # Compute the denominator of the Kazhdan-Lusztig formula
        rho = self.rho
        try:
            # self.load_W_denominator(order)
            self.T_denonminator = self.get_translations_by_n_shift(rho, order)
            self.W_denonminator = self.get_Weyl_group_for_q_series(T=self.T_denonminator)
        except:
            print("\t>>> File not found; creating W_denonminator.", flush=True)
            self.T_denonminator = self.get_translations_by_n_shift(rho, order)
            self.W_denonminator = self.get_Weyl_group_for_q_series(T=self.T_denonminator)
        den = [{w.action(rho) - rho: (-1) ** (w.length() % 2)} for w in self.W_denonminator]
        self.den = den
        print(">>> Denominator created.\n", flush=True)
        return den

    Kazhdan_Lusztig_denominator = CharacterDen

    def Kazhdan_Lusztig(self, order=2):
        print("Computing Kazhdan-Lusztig to order ", order)
        start = time.time()
        q = var("q")

        numerator = self.progress_bar(
            self.num,
            lambda chunk: sum(
                [
                    SR(list(entry.values())[0]).subs({q: 1})
                    * self.character_contribution_from_weight(list(entry.keys())[0])
                    for entry in chunk
                ]
            ),
            100,
        )
        end1 = time.time()
        print(end1 - start)

        denominator = self.progress_bar(
            self.den,
            lambda chunk: sum(
                [
                    SR(list(entry.values())[0]).subs({q: 1})
                    * self.character_contribution_from_weight(list(entry.keys())[0])
                    for entry in chunk
                ]
            ),
            100,
        )

        end2 = time.time()
        print(end2 - end1, flush=True)
        ind = simplify((numerator / denominator).taylor(q, 0, order))
        end = time.time()
        print(end - end2, flush=True)
        print(">>> Completed: %s s" % str(end - start), flush=True)
        return ind

    def SaveNum(self):
        # after calling self.Kazhdan_Lusztig_numerator, save the numerator info into a txt file
        print(">>> Saving numerator to file.")
        num_str = str(self.num).replace(":", ",")
        num_str = "{" + num_str[1:]
        num_str = num_str[:-1] + "}"

        HWString = str(self.llambda).replace("Lambda", "omega").replace("*", "")
        filename = (
            "numerators/num-" + str(self.cartanType).replace(" ", "") + "-(" + HWString + ").txt"
        )
        file = open(filename, "w")
        file.write(num_str)

        file.close()
        self.QSave()
        print("file name = ", filename, "\n")
        return None

    save_numerator = SaveNum

    def SaveDen(self):
        print(">>> Saving denominator to file.")
        # after calling self.Kazhdan_Lusztig_denominator, save the numerator info into a txt file
        den_str = str(self.den).replace(":", ",")
        den_str = "{" + den_str[1:]
        den_str = den_str[:-1] + "}"
        filename = "denominators/den-" + str(self.cartanType).replace(" ", "") + ".txt"
        file2 = open(filename, "w")
        file2.write(den_str)
        file2.close()
        print(">>> file name = ", filename, "\n")
        return None

    save_denominator = SaveDen

    def character_contribution_from_weight(self, weight):
        L = self.weight_lattice
        z = [var("b" + str(i)) for i in range(0, self.r + 1)]
        q = var("q")
        return prod(
            [
                z[i] ** (self.AffineScalarProduct(weight, L.simple_roots()[i]))
                for i in range(1, self.r + 1)
            ]
        ) * q ** (-self.Tolambdakn(weight)[-1])

    def CharContribFromWeightSeries(self, V, weight, w, order):
        # get the contribution from the weights w.(weight - nδ), n = 0, 1, ..., order:
        # - V is the integrable module in question
        # - w is an element of the finite Weyl group
        # A weight with finite Dynkin label [λ1, ..., λr] contributes z1^λ1
        # z2^λ2 ...

        # work out the multiplicity of the weights w, w-δ, w-2δ, ..., w-order*δ
        string = V.strings(order)[weight]
        delta = self.delta
        L = self.weight_lattice
        z = [var("z" + str(i)) for i in range(0, self.rank + 1)]
        max_n = min(order, len(string) - 1)
        return sum(
            [
                string[n]
                * prod(
                    [
                        z[i]
                        ** (
                            self.AffineScalarProduct(
                                w.action(weight - n * delta), L.simple_roots()[i]
                            )
                        )
                        for i in range(1, self.rank)
                    ]
                )
                * q ** (-self.Tolambdakn(w.action(weight - n * delta))[-1])
                for n in range(0, max_n + 1)
            ]
        )

    character_contribution_from_weight_series = CharContribFromWeightSeries

    def CharacterOfIntegrableModule(self, V, order):
        # Return coefficients accurate up to q^order using exact prefix filtering:
        # for each candidate representative, keep it only when its leading q-degree
        # d0 <= order, and then only sum string depth n <= order - d0.
        #
        # Note on Sage strings(depth): the prefix can change when depth increases,
        # so we stabilize only the needed string-prefix (per dominant maximal
        # weight) before assembling the character once.
        q = var("q")
        target_order = int(order)
        if target_order < 0:
            return 0

        wMaxDoms = V.dominant_maximal_weights()
        self.T = list(
            set(
                itertools.chain(
                    *[
                        self.GetTranslationsBynShift(wMaxDom, target_order, order_min=0)
                        for wMaxDom in wMaxDoms
                    ]
                )
            )
        )
        self.W = self.GetWeylGroup(0, "qSeries")

        # Build orbit representatives and compute their leading q-degree d0.
        # d0 = -grade(w.action(weight)).
        reps_by_weight = {}
        max_needed_depth_by_weight = {}
        for weight in wMaxDoms:
            representatives = {}
            for w in self.W:
                acted_weight = w.action(weight)
                key = tuple(acted_weight.to_vector())
                current = representatives.get(key)
                if current is None or w.length() < current.length():
                    representatives[key] = w

            selected_reps = []
            max_needed = -1
            for wrep in representatives.values():
                d0 = -self.Tolambdakn(wrep.action(weight))[-1]
                if d0 > target_order:
                    continue
                nmax = target_order - d0
                if nmax >= 0:
                    selected_reps.append((wrep, int(nmax)))
                    if nmax > max_needed:
                        max_needed = int(nmax)

            reps_by_weight[weight] = selected_reps
            max_needed_depth_by_weight[weight] = max_needed

        global_needed_n = (
            max(max_needed_depth_by_weight.values()) if max_needed_depth_by_weight else -1
        )
        if global_needed_n < 0:
            return 0

        strings_data = self._GetStableStringsPrefix(V, max_needed_depth_by_weight)

        delta = self.delta
        L = self.weight_lattice
        z = [var("z" + str(i)) for i in range(0, self.rank + 1)]

        total = 0
        for weight in wMaxDoms:
            if max_needed_depth_by_weight[weight] < 0:
                continue
            string = strings_data[weight]
            for wrep, nmax in reps_by_weight[weight]:
                upper = min(nmax, len(string) - 1)
                if upper < 0:
                    continue
                total = total + sum(
                    [
                        string[n]
                        * prod(
                            [
                                z[i]
                                ** (
                                    self.AffineScalarProduct(
                                        wrep.action(weight - n * delta), L.simple_roots()[i]
                                    )
                                )
                                for i in range(1, self.rank)
                            ]
                        )
                        * q ** (-self.Tolambdakn(wrep.action(weight - n * delta))[-1])
                        for n in range(0, upper + 1)
                    ]
                )

        return sum([total.coefficient(q, n) * q**n for n in range(0, target_order + 1)])

    def _GetStableStringsPrefix(self, V, required_n_by_weight, max_rounds=8):
        # required_n_by_weight[weight] = max n needed for this dominant maximal weight
        # We only care about those prefixes; increase depth until all required prefixes
        # stop changing.
        required = {w: int(n) for w, n in required_n_by_weight.items() if int(n) >= 0}
        if not required:
            return {}

        # Sage's strings(depth=d) returns d coefficients indexed from 0 to d-1.
        depth = max(required.values()) + 1
        previous_prefix = None
        last_data = None

        for _ in range(max_rounds):
            data = V.strings(depth)
            prefixes = {}
            complete = True
            for w, nmax in required.items():
                if w not in data:
                    complete = False
                    break
                seq = data[w]
                if len(seq) < nmax + 1:
                    complete = False
                    break
                prefixes[w] = tuple(seq[0 : nmax + 1])

            if complete and previous_prefix is not None and prefixes == previous_prefix:
                return data

            if complete:
                previous_prefix = prefixes
                last_data = data

            depth = depth * 2

        if last_data is not None:
            print(
                "Warning: strings prefix did not stabilize within max_rounds; "
                "using latest available prefix.",
                flush=True,
            )
            return last_data

        # Fallback: return the latest call even if incomplete (best effort)
        return V.strings(depth)

    def _CharacterOfIntegrableModuleRaw(self, V, order):
        # internal raw summation for a fixed cutoff order
        wMaxDoms = V.dominant_maximal_weights()
        self.T = list(
            set(
                itertools.chain(
                    *[
                        self.GetTranslationsBynShift(wMaxDom, order)
                        for wMaxDom in V.dominant_maximal_weights()
                    ]
                )
            )
        )
        self.W = self.GetWeylGroup(0, "qSeries")

        # For singular maximal weights, summing over all finite Weyl elements
        # may repeatedly hit the same image weight due to non-trivial stabilizer.
        # Deduplicate by acted weight (within each dominant maximal weight) and
        # keep one representative Weyl element per orbit image.
        total = 0
        for weight in wMaxDoms:
            representatives = {}
            for w in self.W:
                acted_weight = w.action(weight)
                key = tuple(acted_weight.to_vector())
                current = representatives.get(key)
                if current is None or w.length() < current.length():
                    representatives[key] = w
            total = total + sum(
                [
                    self.CharContribFromWeightSeries(V, weight, wrep, order)
                    for wrep in representatives.values()
                ]
            )

        return total

    # alias
    character_of_integrable_module = CharacterOfIntegrableModule

    def progress_bar(self, full_list, chunk_processor, chunk_size):
        # perform the first computation
        # use to initialize the variable <result>, which can have different
        # data structure depending on the task
        result = None
        for i in tqdm(range(0, len(full_list), chunk_size), leave=True):
            chunk = full_list[i : i + chunk_size]
            if result is None:
                result = chunk_processor(chunk)
            else:
                result = result + chunk_processor(chunk)
        return result
