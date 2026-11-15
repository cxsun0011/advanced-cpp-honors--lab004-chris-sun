#include "project/DynamicProgramming.hpp"
#include <algorithm>
#include <climits>
#include <stdexcept>

std::string buildParenthesization(
    const std::vector<std::vector<int>>& split,
    int i,
    int j) {

    // TODO: Implement recursive reconstruction of the optimal
    // matrix-chain parenthesization using the split table.
    // Base case: if i == j, return "A" + std::to_string(i).
    // Recursive case: split at split[i][j] and combine both subchains.

   if (i == j) {
        return "A" + std::to_string(i);
    }
    
    int k = split[i][j];
    std::string leftStr = buildParenthesization(split, i, k);
    std::string rightStr = buildParenthesization(split, k + 1, j);
    
    return "(" + leftStr + " x " + rightStr + ")";
}

MatrixChainResult matrixChainOrder(
    const std::vector<int>& dimensions) {

    // TODO: Implement Matrix Chain-Product using bottom-up dynamic programming.
    // dimensions = {p0, p1, ..., pn} represents n matrices.
    // Create a cost table and a split table, evaluate increasing chain lengths,
    // store the best split, and reconstruct the optimal parenthesization.
    // Validate that at least one matrix exists and all dimensions are positive.
    // Do not hard-code the sample answer.

    if (dimensions.empty() || dimensions.size() < 2) {
        throw std::invalid_argument("Invalid dimensions vector for matrix chain multiplication");
    }

    for (int dim : dimensions) {
        if (dim <= 0) {
            throw std::invalid_argument("Matrix dimensions must be positive integers");
        }
    }

    int n = static_cast<int>(dimensions.size()) - 1; // Number of matrices

    // cost[i][j] stores the minimum scalar multiplications needed to compute matrix A[i..j]
    std::vector<std::vector<long long>> cost(n + 1, std::vector<long long>(n + 1, 0));
    
    // split[i][j] stores the optimal index k where the subchain is split
    std::vector<std::vector<int>> split(n + 1, std::vector<int>(n + 1, 0));

    // l is the chain length
    for (int l = 2; l <= n; ++l) {
        for (int i = 1; i <= n - l + 1; ++i) {
            int j = i + l - 1;
            cost[i][j] = LLONG_MAX;

            for (int k = i; k <= j - 1; ++k) {
                // cost = cost of multiplying left subchain + cost of right subchain + multiplication cost of the two resulting matrices
                long long q = cost[i][k] + cost[k + 1][j] + static_cast<long long>(dimensions[i - 1]) * dimensions[k] * dimensions[j];
                
                if (q < cost[i][j]) {
                    cost[i][j] = q;
                    split[i][j] = k;
                }
            }
        }
    }

    std::string paren = buildParenthesization(split, 1, n);

    return MatrixChainResult{
        cost[1][n],
        split,
        paren
    };
}

AlignmentResult globalSequenceAlignment(
    const std::string& first,
    const std::string& second,
    int matchScore,
    int mismatchScore,
    int gapScore) {

    // TODO: Implement global sequence alignment using dynamic programming.
    // Build and initialize the DP table, evaluate diagonal/up/left transitions,
    // then backtrack from the final cell to reconstruct both aligned strings.
    // Use '-' for gaps and support arbitrary scoring parameters.
    // Do not hard-code the DNA/text examples.

    int m = static_cast<int>(first.length());
    int n = static_cast<int>(second.length());

    // dp[i][j] stores the maximum alignment score for first[0..i-1] and second[0..j-1]
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));

    // Initialize base cases (gaps along the edges)
    for (int i = 0; i <= m; ++i) {
        dp[i][0] = i * gapScore;
    }
    for (int j = 0; j <= n; ++j) {
        dp[0][j] = j * gapScore;
    }

    // Fill DP table
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int scoreDiag = dp[i - 1][j - 1] + (first[i - 1] == second[j - 1] ? matchScore : mismatchScore);
            int scoreUp = dp[i - 1][j] + gapScore;
            int scoreLeft = dp[i][j - 1] + gapScore;

            dp[i][j] = std::max({scoreDiag, scoreUp, scoreLeft});
        }
    }

    // Backtrack to reconstruct the aligned sequences
    std::string alignedFirst = "";
    std::string alignedSecond = "";
    int i = m;
    int j = n;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            int matchOrMismatch = (first[i - 1] == second[j - 1] ? matchScore : mismatchScore);
            if (dp[i][j] == dp[i - 1][j - 1] + matchOrMismatch) {
                alignedFirst += first[i - 1];
                alignedSecond += second[j - 1];
                --i;
                --j;
                continue;
            }
        }

        if (i > 0 && dp[i][j] == dp[i - 1][j] + gapScore) {
            alignedFirst += first[i - 1];
            alignedSecond += '-';
            --i;
        } else {
            alignedFirst += '-';
            alignedSecond += second[j - 1];
            --j;
        }
    }

    std::reverse(alignedFirst.begin(), alignedFirst.end());
    std::reverse(alignedSecond.begin(), alignedSecond.end());

    return AlignmentResult{
        dp[m][n],
        alignedFirst,
        alignedSecond
    };
}
