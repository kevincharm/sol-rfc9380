// SPDX-License-Identifier: MIT
pragma solidity ^0.8;

import {ModexpInverse, ModexpSqrt} from "./ModExp.sol";

/// @title SVDW
/// @notice Shallue-van de Woestijne mapping for BN254 (G1)
contract SVDW {
    // Field order
    uint256 private constant N =
        21888242871839275222246405745257275088696311157297823662689037894645226208583;

    /// @notice Param A of BN254
    uint256 private constant A = 0;
    /// @notice Param B of BN254
    uint256 private constant B = 3;
    /// @notice Param Z for SVDW over E
    uint256 private constant Z = 1;
    /// @notice g(Z) where g(x) = x^3 + 3
    uint256 private constant C1 = 0x4;
    /// @notice -Z / 2 (mod N)
    uint256 private constant C2 =
        0x183227397098d014dc2822db40c0ac2ecbc0b548b438e5469e10460b6c3e7ea3;
    /// @notice C3 = sqrt(-g(Z) * (3 * Z^2 + 4 * A)) (mod N)
    ///     and sgn0(C3) == 0
    uint256 private constant C3 =
        0x16789af3a83522eb353c98fc6b36d713d5d8d1cc5dffffffa;
    /// @notice 4 * -g(Z) / (3 * Z^2 + 4 * A) (mod N)
    uint256 private constant C4 =
        0x10216f7ba065e00de81ac1e7808072c9dd2b2385cd7b438469602eb24829a9bd;

    error InvalidFieldElement(uint256 x);
    error MapToPointFailed(uint256 noSqrt);

    /// @notice Map field element to E using SvdW
    /// @param u Field element to map
    /// @return p Point on curve
    function mapToPoint(uint256 u) public view returns (uint256[2] memory p) {
        if (u >= N) revert InvalidFieldElement(u);

        uint256 tv1 = mulmod(mulmod(u, u, N), C1, N);
        uint256 tv2 = addmod(1, tv1, N);
        tv1 = addmod(1, N - tv1, N);
        uint256 tv3 = inverse(mulmod(tv1, tv2, N));
        uint256 tv5 = mulmod(mulmod(mulmod(u, tv1, N), tv3, N), C3, N);
        uint256 x1 = addmod(C2, N - tv5, N);
        uint256 x2 = addmod(C2, tv5, N);
        uint256 tv7 = mulmod(tv2, tv2, N);
        uint256 tv8 = mulmod(tv7, tv3, N);
        uint256 x3 = addmod(Z, mulmod(C4, mulmod(tv8, tv8, N), N), N);

        bool hasRoot;
        uint256 gx;
        if (legendre(g(x1)) == 1) {
            p[0] = x1;
            gx = g(x1);
            (p[1], hasRoot) = sqrt(gx);
            if (!hasRoot) revert MapToPointFailed(gx);
        } else if (legendre(g(x2)) == 1) {
            p[0] = x2;
            gx = g(x2);
            (p[1], hasRoot) = sqrt(gx);
            if (!hasRoot) revert MapToPointFailed(gx);
        } else {
            p[0] = x3;
            gx = g(x3);
            (p[1], hasRoot) = sqrt(gx);
            if (!hasRoot) revert MapToPointFailed(gx);
        }
        if (sgn0(u) != sgn0(p[1])) {
            p[1] = N - p[1];
        }
    }

    /// @notice g(x) = y^2 = x^3 + 3
    function g(uint256 x) private pure returns (uint256) {
        return addmod(mulmod(mulmod(x, x, N), x, N), B, N);
    }

    function cmov(uint256 x, uint256 y, bool b) private pure returns (uint256) {
        if (b) {
            return y;
        }
        return x;
    }

    function sgn0(uint256 x) private pure returns (uint256) {
        return x % 2;
    }

    /// @notice Compute Legendre symbol of u
    /// @param u Field element
    /// @return 1 if u is a quadratic residue, -1 if not, or 0 if u = 0 (mod p)
    function legendre(uint256 u) private view returns (int8) {
        (bool success, bytes memory _x) = modexpLegendre(u);
        if (!success) revert MapToPointFailed(u);
        uint256 x = uint256(bytes32(_x));
        if (x == N - 1) {
            return -1;
        }
        if (x != 0 && x != 1) {
            revert MapToPointFailed(u);
        }
        return int8(int256(x));
    }

    /// @notice This is cheaper than an addchain for exponent (N-1)/2
    function modexpLegendre(
        uint256 u
    ) private view returns (bool success, bytes memory output) {
        output = new bytes(32);
        bytes memory input = new bytes(192);
        assembly {
            let p := add(input, 32)
            mstore(p, 32) // len(u)
            p := add(p, 32)
            mstore(p, 32) // len(exp)
            p := add(p, 32)
            mstore(p, 32) // len(mod)
            p := add(p, 32)
            mstore(p, u) // u
            p := add(p, 32)
            mstore(
                p,
                0x183227397098d014dc2822db40c0ac2ecbc0b548b438e5469e10460b6c3e7ea3
            ) // (N-1)/2
            p := add(p, 32)
            mstore(
                p,
                0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
            ) // N

            success := staticcall(
                gas(),
                5,
                add(input, 32),
                192,
                add(output, 32),
                32
            )
        }
    }

    /// @notice sqrt(xx) mod N
    /// @param xx Input
    function sqrt(uint256 xx) internal pure returns (uint256 x, bool hasRoot) {
        x = ModexpSqrt.run(xx);
        hasRoot = mulmod(x, x, N) == xx;
    }

    /// @notice a^{-1} mod N
    /// @param a Input
    function inverse(uint256 a) internal pure returns (uint256) {
        return ModexpInverse.run(a);
    }
}
