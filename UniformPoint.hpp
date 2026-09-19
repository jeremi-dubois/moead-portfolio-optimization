#pragma once
#include <vector>
#include <cstdlib>
#include <utility>

// ============================================================================
// Génération de points uniformes sur le simplexe unité (méthode Das-Dennis /
// simplex-lattice), pour initialiser les vecteurs de poids/référence.
//
// H est choisi automatiquement pour que le nombre de points combinatoires
// C(H+M-1, M-1) soit le plus proche possible de N demandé.
// ============================================================================


// Calcule le coefficient binomiale C(n, k)
inline long long binomCoeff(int n, int k) {
    if (k < 0 || k > n) return 0;
    long long res = 1;
    for (int i = 0; i < k; ++i) {
        res *= (n - i) / (i + 1);
    }
    return res;
}


// Génère les combinaisons récursivement (c_1,..., c_M) telles que sum(c_i) = H, c_i >= 0, 
// puis normalise par H pour obtenir des poids sur le simplexe
// M         : nb objectifs (M > 0)
// H         : nb subdivision du simplexe (H > 0)
// remaining : somme restante à distribuer entre les composantes
// idx       : indice de la composante actuellement générée
// current   : combinaison entière en cours de construction
// out       : ensemble des vecteurs de poids générés
inline void generateSimplexLatticeRec(int M, int H, int remaining, int idx,
                                       std::vector<int>& current,
                                       std::vector<std::vector<double>>& out) {
    // Cas de base : dernière composante atteinte.
    if (idx == M - 1) {
        current[idx] = remaining;
        std::vector<double> point(M);
        for (int i = 0; i < M; ++i) point[i] = static_cast<double>(current[i]) / H;
        out.push_back(point);
        return;
    }

    for (int c = 0; c <= remaining; ++c) {
        current[idx] = c;
        generateSimplexLatticeRec(M, H, remaining - c, idx + 1, current, out);
    }
}


inline std::pair<std::vector<std::vector<double>>, int> uniformPoint(int N, int M) {
    // On a donc deux nombres de points qui encadrent N : C(H+M-1, M-1) <= N et C(H+M, M-1) > N
    int H = 1;
    while (binomCoeff(H + M, M - 1) <= N) H++;

    // Compare H et H-1 pour choisir le plus proche
    long long countH   = binomCoeff(H + M - 1, M - 1);
    long long countHm1 = binomCoeff(H - 1 + M - 1, M - 1);
    if (llabs(countHm1 - N) <= llabs(countH - N)) H = H - 1;

    if (H < 1) H = 1; // Sécurité

    std::vector<std::vector<double>> W;
    std::vector<int> current(M, 0);
    generateSimplexLatticeRec(M, H, H, 0, current, W);

    return {W, static_cast<int>(W.size())};
}
