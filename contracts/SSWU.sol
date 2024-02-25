// SPDX-License-Identifier: MIT
pragma solidity ^0.8;

import {ModexpInverse, ModexpSqrt} from "./ModExp.sol";

/// @title SSWU
/// @notice Simplified Shallue-van de Woestijne-Ulas mapping for BN254 (G1)
contract SSWU {
    // Field order
    // prettier-ignore
    uint256 private constant N = 21888242871839275222246405745257275088696311157297823662689037894645226208583;

    /// @notice Param A of E', a curve isogenous to BN254
    uint256 private constant A =
        9087994317191712533568698403530528306233527979934880849865820425505218365052;
    /// @notice Param B of E', a curve isogenous to BN254
    uint256 private constant B =
        3059101143800926337153883959975852125336293569895750485959800095292563537400;
    /// @notice Param Z for SSWU over E'
    uint256 private constant Z =
        21888242871839275222246405745257275088696311157297823662689037894645226208570;
    /// @notice (N - 3) / 4
    uint256 private constant C1 =
        5472060717959818805561601436314318772174077789324455915672259473661306552145;
    /// @notice sqrt(-Z)
    uint256 private constant C2 =
        9128901577248916245119121441495541088906573436435020319814006418557303224083;

    error InvalidFieldElement(uint256 x);
    error MapToPointFailed(uint256 noSqrt);

    /// @notice Map field element to E using SSWU
    /// @param u Field element to map
    /// @return p Point on curve
    function mapToPoint(uint256 u) public view returns (uint256[2] memory p) {
        if (u >= N) revert InvalidFieldElement(u);

        uint256 tv1 = mulmod(u, u, N);
        tv1 = mulmod(Z, tv1, N);
        uint256 tv2 = mulmod(tv1, tv1, N);
        tv2 = addmod(tv2, tv1, N);
        uint256 tv3 = addmod(tv2, 1, N);
        tv3 = mulmod(B, tv3, N);
        uint256 tv4 = cmov(Z, N - tv2, tv2 != 0);
        tv4 = mulmod(A, tv4, N);
        tv2 = mulmod(tv3, tv3, N);
        uint256 tv6 = mulmod(tv4, tv4, N);
        uint256 tv5 = mulmod(A, tv6, N);
        tv2 = addmod(tv2, tv5, N);
        tv2 = mulmod(tv2, tv3, N);
        tv6 = mulmod(tv6, tv4, N);
        tv5 = mulmod(B, tv6, N);
        tv2 = addmod(tv2, tv5, N);
        uint256 x = mulmod(tv1, tv3, N);
        (bool is_gx1_square, uint256 y1) = sqrtRatio(tv2, tv6);
        uint256 y = mulmod(tv1, u, N);
        y = mulmod(y, y1, N);
        x = cmov(x, tv3, is_gx1_square);
        y = cmov(y, y1, is_gx1_square);
        bool e1 = sgn0(u) == sgn0(y);
        p[1] = cmov(N - y, y, e1);
        uint256 tv4_inv = inverse(tv4);
        p[0] = mulmod(x, tv4_inv, N);

        (p[0], p[1]) = isoMap(p[0], p[1]);
    }

    /// @notice TODO: Replace with addition chain
    function modexp(
        bytes memory base,
        bytes memory exponent,
        bytes memory modulus
    ) private view returns (bool success, bytes memory output) {
        bytes memory input = abi.encodePacked(
            uint256(base.length),
            uint256(exponent.length),
            uint256(modulus.length),
            base,
            exponent,
            modulus
        );

        output = new bytes(modulus.length);

        assembly {
            success := staticcall(
                gas(),
                5,
                add(input, 32),
                mload(input),
                add(output, 32),
                mload(modulus)
            )
        }
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

    /// @notice Optimised sqrt_ratio for curve order equivalent to 3 (mod 4)
    function sqrtRatio(
        uint256 u,
        uint256 v
    ) private view returns (bool isQR, uint256 y) {
        uint256 tv1 = mulmod(v, v, N);
        uint256 tv2 = mulmod(u, v, N);
        tv1 = mulmod(tv1, tv2, N);
        // TODO: Exponent ^C1 should be precomputed
        // This currently takes about ~2450 gas
        (bool succ, bytes memory y1_bytes) = modexp(
            abi.encodePacked(tv1),
            abi.encodePacked(C1),
            abi.encodePacked(N)
        );
        require(succ, "modexp failed");
        uint256 y1 = uint256(bytes32(y1_bytes));
        y1 = mulmod(y1, tv2, N);
        uint256 y2 = mulmod(y1, C2, N);
        uint256 tv3 = mulmod(y1, y1, N);
        tv3 = mulmod(tv3, v, N);
        isQR = tv3 == u;
        y = cmov(y2, y1, isQR);

        return (isQR, y);
    }

    /// @notice a^{-1} mod N
    /// @param a Input
    function inverse(uint256 a) internal pure returns (uint256) {
        return ModexpInverse.run(a);
    }

    /// @notice 59-isogeny map E'(x',y') -> E(x,y)
    /// @param xPrime x'
    /// @param yPrime y'
    function isoMap(
        uint256 xPrime,
        uint256 yPrime
    ) internal pure returns (uint256 x, uint256 y) {
        uint256 x_num;
        uint256 y_num;
        uint256 x_den;
        uint256 y_den;
        assembly {
            let x_i := 1
            x_num := addmod(
                x_num,
                mulmod(
                    0x12657257b36598525b8dedc31e3dcbc917f6242bd32bf287d349a3075a9ae6c4,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x3048b9d6aee10c2ae3f7f100bf3595069a5edb5aa8079f3ba4a365d9b334da82,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0xba58d92e0b040dc8fc21576230e7f21f85025d2e75f37745f6b16a3701f458,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1cf56f40dd3bb57fd5108c04253f511630e3ec8af7eec55a6de515f8ee0395aa,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x22571421f593331f844eabf61502d36a169475ff48973acf17485172b67b32d2,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x2ec29c2b03452b83e591f48d047e2342f203129817fc82d2753e37e8dc38564,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x2b419ea24f5be6297ea7ba4def319adace88089a99554cc2f7f280cbef3987f0,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1e6d2100e2a9d281e1ce627f18fdaac39a275f31e7539fe477b7ca6a9fabc791,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0xa512f2f0f576d2454332062a8a2b00f46145c24dda3e0328f0173c38c88f28c,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x120a4376efb9d2da3e6465e9ff32fa0ad96632dfbb1b29cce48e9d5a30e73d41,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x2d9dd84348d9c18d045c4c5570e69426934669180054657939ab3a302b9d6758,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x23dd6c88f9ea285877b582475685cd45928a04b312d5b2e0b096c3784ef70948,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x2236c6aaa098e8057f81a8f490e5c760664d53b540e7074c1763b459ced63fa8,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x63949278e29d308dda0aaffd95962b0310fa6d6a28259b9d5f19c20ce3ce20b,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1a8928be12dc6d0534053b5f48ce308455b384f73724c8f15e43438a38303acc,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2af5ba4ca63607822516c7e3350cebff03b100234cf91321f502074807ac8a52,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x248dadf66b33a1c97393709d3548c9abac08e9941a76c998b41b992b2c2e0f2e,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1b7cc3a7950f4fb21b381e1f24d1114fa225c5883cf907d45b6a4793c91a40d0,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x293c397ef7065152e3345dacc78e7953230197f93cb441d0e47d50c6a0597830,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x5b71c393b7b4c7fd9b16e59b95d1de70d87614ec9c07796ed089c53a353eb7,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x679488d8950a9e8988cd9b3e7a2b25c511fdf49f47c0de308d2083d5766cc7,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x134ca60f6e1a628c013ffbebd0ef1e005cdffcd039d7e2fc1435a85073178f1a,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x142099a95adb4dea32e513b0e17bbb25adcb28e3ae6b47827ebc6145d49d3e7,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x10525931ea4f68c977503426d087ce7a3bec5915ded132f06ffe713d2b2a49e0,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x14e939bac3ce8bfb846632f83393dbbb58ad309e63d8b86483256bc0b6bf1c81,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x431f8a7b99b00318c3b82c67ec2471c0382f5e53ce91c4c4893467821274331,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0xa954e212071f50714810ac25fe967cf850e0f3329ef104ab15ad1881b19b050,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2848a31144faeb433aa793e8f26acff1b523db795880ba3ce59756daee27833e,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x9204710f6347df5b7a2b467f2eb4033215d4dce3a7f51d3c6502f345e4e851a,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0xead42c16d9428d995e99602232d467bd1286df3a0dac44f9f0696603647d8ae,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x8b44c6e1990189af55373c9acf347369b7abad5a131e87bccc0cbcd077396dd,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0xa3f7029cfe58bfd3e5cc8369a39c90b04d25f567bc174bc12234466358bb0d5,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x24295418cb5ff8285aa3f5fe44a4c24f936c1016e1b44581c81c7bc45947b0d,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x19fe81719193f268ef87815f78237ad2b0815fdb8002dc59e63234765d9fff93,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x19ccbf6a84083d9f9cb43338258082e08a6e78da7d147f7a52a62a4b6418d187,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x3dd73ed02aa13f476d41eabf74197fc49389739f1bdb9df9b323515e8dad6f9,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1236daf5fa1eab478b090b4441d70698a7512b5ad22a25253b6165a9dc82f93,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x240b333d2efafd3ebacb6f87338bc18d7c03d74838642de21398d71ff39e877d,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x123578aca9c6ba8c5fe428a505f0f1ba324a144f2ed31fb64031b3b0cbe3f5d5,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1ae75965a5d2d4319b6bb44a689a5093b3d5a84c0ab6bceef529870495cda742,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x122d907f66040c6cd6ee369c90ff243a5661f65ce4a91015aa01a160352788a3,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0xa220b3ad98bb608bd64fe0595fe6e906ef1b792764fb47a258be0efa17f961e,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x5fa3045ff01f64dec739aa0cf87d9f960e45223dc753508de3cc9f53c3d81d6,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x14664f4ca095030dab68dc49a4213ac8ae13c2e76f8572f0d1a014f09ad46e83,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0xcbaaef937f31477029bcfa1237757e2c153e463b15d58fa5927c0dee0c0455c,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x23bab3ea1049c6c97a236e5cd423c91a7238a68d50bf401746831ccf1346fb0b,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x2dcf52a50488060b325ccff059f2899de781071560632626462ff5599a107196,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2789f9d7077057210a526cd335794b4140e31809a756990dafabc85b998b6dc5,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x2d605fba4020d945c2e74278890bda29643bd021a61b4c3deeb5b195126aa2a1,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x8247b9e325a2ac12fe815df3dfc6cb962d68f82aaabe0fb20eed0ef99e9e172,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x2b093063280bc9a335b9283ccd6db69ee80b23f297b3b0bc0dc9d4dc3c56f002,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x5cf4627f6525317a6456e150531e566736a7bae4c4a550c97270f8e668a8cd4,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x131bd55912887e270454b5bf07fe816ac3753e16b07f712115720b8502363532,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x11c1d2dc85853afdda57a325f0efbf3a79dacc4eef29fd168b0ad9bb65836cfe,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0xaa48ae9dd0833514a969c7791afbe4a2670532caa14b6a094d4825064d7a37f,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x380855421858e14bf09a17b6279985200c942d8be15ab34e26c66123ac2afba,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0xa5b43bf861477b9677f327adfff4b540e3beaeed33acf9ece454441908c416,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1d35df5e1ae75b4f98cdf390e6d154d21a9e7bc0fe0b0db4e42b757f1c49be27,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x29e5ac563f3de053539bff1a2f687c622574dd49abe1973371ca050c78193fde,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x25b81ce31f172aa593c7eb3a9063af5bd2a4a39c8cdb597a373f9f6991117c37,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x43756a8089a53919be72e96ae9a2d0d7d9a7b78bd2ac4f08be4e49f3f2496aa,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x17502514fd67654bc4e19144e8db83492652c2ec4f1df4f06ea23aef90949d25,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1d3ff29dd0345d9250457fc26222f1b2cb5140329b970b9ce06cd001786eb3a9,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1ac06996e9de305842b6b5982fe25c1359de74e0d2c84dc374aa764c8047d323,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x204c4583f0fdd45b7ddda85cca9b5968a5bbbf2a0fafc202f5b0b6334c335645,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x996c7e6686e305cf4e7dfb1bc3ca5d54125eca241996d01fcc4def2e0f4c864,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x13d100a74a8b38cc399eb8a23362acb3fb4a6ae15921e204a6fa2e125852c3f,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0xdcf3094d55e0f87900e1cab28dd7d9cad69d7d9fc1e8ae8c62841a423269d34,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x163af0e1c92c8412575e76ac6e96c55998e8c9f3710127a22b12aa56fd77c3e1,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x178907477e537a01f6f40a0dc80efadc10c5bd51624ae31bf55ebeb71df15a16,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x5d60597ff6f7cfd31cfa568d4dad01f50eaca7c79d7e2b62c8d91ff63092ee8,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2a7e4f63cf5531bb934388c8be0797a634d0df9064efceb8a3772c6455f3cdfd,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1e3e80d7d73787f723e228e9a23127a5b7c1ecf59a34945c2cbfffe87efc96c7,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x2a61c00b278ecc294a73cc7c6febf717775dbee86b33e35ba272cdca86067e0d,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x22a9ed5289c350ae8d9db43eaa80eab52eaf4c1ac5db3f66cd42c8de6eabaff0,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x763f3fe0cd4bab574d5f7d91edf5ccde03ff7b6d1ee88cb6a26e2c58e5ded7e,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1aa0c6e2b90c3fb8a2907372eb0bca814b0d4eea342dfc093f382ca0a28f4623,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1b03b36b77c7c28e8344684c650e54225f770980b237e93754b1179b030741ee,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x2997a6c02b67c9f32aa1e4ead96574cb1c94b9846489419552ce45e7759e0235,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x159df2d2b09ac96cf08bd103d6cc55ec148111fc6e66d089a19ce00bb5d4555b,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x6a781fa1db98fa7662a38112455e5a6364c76efdcfee81c02394c3944715905,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1f0d4f7e131580b0b4057705f971179831cccb30ad9c98808c0cc898882f1cce,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1a139a3da3c49bb1a4314ca87e8a99f7861dc840d11476cdf12dd2887f9da6c6,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x232504197d88a44c849fcd648f8974bfd488ff51d77a4a57f803fff1648d0362,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x27e79fdf14e31c92255208fe43e328a100c29c3b169728dae753913037259f29,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0xb0b56ebe38e2ae08729838fb12b537da910feedbec6d3f2fb56b797e7a535f7,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x7f2077d6d88b3798673ba35cf18b9b77dce8a19e823df93a98aae54c9a20f72,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x14f06c9d480abcef29196bf0e02ce5723475e6fa4ae4fdc4b5d1036a65f7eb14,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x167c6d2aabe07bced4d43299c7471759d23a2e7fb03829b552af75e2b84cb9a,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x101e21fc2aaf60a9b968612cfd698724cdddae783b9298b248dcbf07b8b37e1f,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x14886addd7135a87b3ed11dac44571f809ca55f778f99e4895ed3f0191b68485,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x10a7d368586bb03c512372d0be40222207c387c05d95d152e1f2486241d059d9,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x2e011bbb4846efa4646c96662f590f2d3391fb048247b9fe6cbbc86ec57422f4,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x18a7e793e1e895cde67246cf55217fcca5715a554cfe80d0de42e117ffdbe0dc,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x248600e850e726404c1369bb80d09dd7a046fba510402bf7bcc4ce2b76b6d887,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x124cf290f9d08ce3813489926386192311afbe119f735c2c4709a2a63581d6b4,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1f37e8b18f561800cb92a7d8b43bae4cdfebcf1f3222bf409f9f44d86362114c,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x7a9029e7beda9d21865830a7e55b99feb30a0860ed5278f2658d2973fff9ff9,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0xc514ba57dd7eb80b4fc3fad1016f13616bff11022534f7ad81178143113b75f,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x24dac56344dbe2bca64911797833f0c95cb3e9d89e985a1b78cb36ef1d3a356c,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x26a2f65fdf1e765eea5c9bc665ce9bd294cdb4231e045b6e8b105d8b205a2cb9,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x21d825e4363b23a8f79130fff89f343a98890fb18f0364ea8b1bfdbb537898d2,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x65fcf803fa1b53e95394ff6ade54648a5fc137b51cf42348420ce30b3bd691e,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x12a3837deaabe36c27e6545067c5e335a5f77b879566d996eb080ac058002815,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x222572c371692ce6782a8fe01af708d995280f6ca355c9270894c6f8cabce4ed,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x8e64c0b02afdb053d9c4e33f44ea2c456c8ba94f18ae5e65a80a54069ed7fdf,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x168dc422ad982df209ea7353065a19d203bca255426cbbb46152081b7d99ff42,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x408569894e16438d8f2aa4d7455131a45e942cac892addbedd5f653f596c688,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0xed07b5e5cf3da8dae7d0201a114bb5a14955bd24b69a77650415a12ec390a79,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0xcfd0c28324dc74fad12eb0ec318ba744abf5531a913839b4f51363d5b0f54e4,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1dad59efda1367c9a23bb12742c0d0f201678ecc123cc88d5994d58c31056c57,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0xcfca66fc31e36b0ac109c6e68626d9ee1823e3ca1198830122a063ae7e1b294,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0xb8adfc72afaafd2c77d1b5d79f0489ffabbdc5c245970432e11b7275a6f2a57,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x49607ed54a606fd365a064985e38bf901bed501a32bc616f71a5f6dca1c69ae,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x188c74b27f5461cfe3425e441965a54b2702d2b0627e4302bdec63e3f80897ee,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2a3fd85c37cee4ad33ab5723d0a464c621df825a3601d12eee4fa247342e55f0,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1656afd683b7ab47562bcd88f31f90faf8bd09ad9004e3fc911d073fb38dc40b,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x971114c1326592948d89e604ee2e31921bea2704f6b02d3bfa86a635f1c5ea5,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x18d75442497967a12066032c7ebee8f6200c7f3507d5a6570c5d1a2ee0c2d631,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2bd9a3129bcd5e6921895fce47ec906050bd9162b57eff29a104f12a77078dab,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x121bab4a9b6b7006d71d334064237e94be4acfcd6ddeddcf47f463b1edbc7be2,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1517e31fa906b7f6714037121b3e992039310302a5e4bd5e09c7030b1355f6d5,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x9356b07ded28a21aa0f2d79f27fcf30c5c5e8e3d16489ada95a40a8e4def4f3,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1cad4c36823efc2445fe46ee20d1dc484d2539a11182d926e296e62cac6f3723,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x9ee0857b3065ea5935bf1ed98c4c1d9eeee73562ee6c7192029c014d59526d7,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1a6f76a7829853e81c150e878be4fe4fc8dba5891d2567230be885baaca0bcf,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1718356a1b0cbadbe436895129a568dde93463b5f3813dee8aaba615161fdfa8,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1dcd7bcc7fcd822b8c34a567ceaa545aba95977014e7c3924261c69824ce8088,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1ad21f50d6bd4201da68e9de35590c54b5178c68685236c024ef1739e6489e9,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x11681c18663a229d284c822d3f2edfff105d90a9cfcbbba3b17d793aa22ecfc7,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1734021f0745afaedd0dc66abb0c91c795114b51846fbf9be9eaf342e5683a4,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x344e387dcb01159881ff3806293311e90762e0c23b6c27a69d4dbed824d4396,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1a810c88f1b8d7f9b9bd64e228ecddd826d6cef79c28ba1c44aecfc54585ba56,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x25e9e7d14f8b6e1dbfe038c50941fcbc8ce2be6f0a485f055abadfb738d8525,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1eadd045a2bd0e2d51ea04fefd335d7f200fc6157c30f04c035fe7037a504b6e,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x27a47369346acb467a6c58b5c4ac5bc720927b36fac7ebe44896d6ce2cb5cf91,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x118492813eecec14d1dfb05e93f7a5ee223369ed9450b4eb3d4b250d2455e344,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x2042e29c879af2171753d85c5131eb38591a11790201c2f00cec4c85c18617b4,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1fbeae0787fe09f6672e16b14711ca04aeab89db424dffbb1891336cd9b74bd4,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x20c04eb5c66677fa0ac834c377e0430ba6a0296512d2851d2a16c39d417b782,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x15be68db677076a8151e4c7236637291d7b82fbca419b95ea5be784f31e1358d,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1ad36c25f49005c8751e9fb60f4bdaf295e199e316fdfafba75ab85c7e8971cf,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1234b827a5f3461cbf0ce2decdfc6beced65002fe8f331c648661e66bf161b0e,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2ca1ed33b3ecd3479279315c236dbb922780d457a0a72a17690f09c8cda0dbe,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1aa3dcb87a76b10973eca16b6b2fd2116069f6ea89c7738158817dafe5b43bd4,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0xbe9c1120726e2f56228f2e2fb11cbba0a426d291d18671f3a65bc8190125270,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x24bde858c13f232450b750253bfa329246c36999407ac581390ffaa2b6acb090,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x9a836b2edad29ca347e8f15c9dc7f5177461a5d8edec335216181ce70ea2cf1,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x60c71d270a49a188f4900239f4b9742f526f7369d7954d9ceb3588b8bedad69,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1e671427d328b62e267a400870840eec56b123c6225859fb21587b3927052f14,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x29026b7b766e4c5c2eeec44b9ff1a5fec10a39888fde90cc702d69f710b44f34,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x8d533e4fdfb3cb1d1db52282ae89a779f6f23d4ae3dd7fd28d4e03860ad13b9,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x2efb45aa203fde80616e09ab55bbbea5cee71812866b932e284bdd21a7e8bb61,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1f7f83cbfcb05b98d570384be701f65db9b8e440d20e5361d2565805361dbcdc,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1a396d53cdb66c65cd2be99442ebdcea553ea7273bed98c48d8f696449aa6306,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1e19abb525eed382cc6eaab21c30289c84ce31c9916b587560c4b337787fab7f,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x2913cdec8bbc1a803a1fdcf954713b9b2f50e2a6eac56392923114b15394e7c7,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1a14e1adda34c1543eb3ed03a47a2973ef9e2bd4da2bdab3c98af414703ecc43,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x98c9574c36d36419869f0f04fb43e93d53dc647a200b5b83302d935c76ef9dd,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1fea36cdb44b8946e8741739d3547bb27ca54be13e61a89c180684d4ac17ce23,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x46c5668cc1ad3a87e573ff3e37726b1688281c64e31d49d4fc12d6559e7f014,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x23420d6be0f292bef0ec34ad3927070c6c9d3b6e6fff8326d896b87ca9397713,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x74b59da1754c36a849895e4e4d160677b8b10d8b152c7cff5f2ae5da667cc3,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x8491b09ae81f55d03c0ffe15af73e3c76316617a8f654e06a18cd1a1e29ba2a,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1e00a0d18c7edd49c3058bbaaba31e82f0ebc1a1c181394a8b1be5768f6e7386,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x26c605120dc11ba4b41d07d57f32efc256dc8b9d8dcf1953f400d2310a476b88,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x15c0b0249f4ec9871bac3db7e7f4ef5e137a1abad24bf8943ffd597c3931e75a,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1383b102615f47d9cb5b835fc18de4ca8826fdc8f745292ed24bde0c2600d0d8,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1b56578e84300fb95174ec6d6969285759c18481b33079d0571bba419c6c23e9,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x282ba6992602793ce135023f205426f86dbb04ee0df6a77a49c17d979a5b6deb,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x7f3d6ee6b19d0bc10bea26cdce5b84d7f4a321d18dd01cafdd56c9e085d1f54,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2939940dfb4b2725d35a35f19b388b99053648cf22822ef9d33b744fbc02e460,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x2648a2a5f5bd2b2fe89121490c99566813f8e6dc58f43f894f2be48fc731fd9f,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1cf40279367084248032d68d8b27dd18b1f17e191b1b9ff0097d053308c185ba,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x27fb26e6f44cc3f5ff9db5b88ce5fc315717b8ccdf7cb59e6ecb1118864c26ad,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x4f283c6f46edd5092562916887a47977dfdf2780d16302e53a2b76ce86d6ed0,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x16df1b875e4c8a5f596953fe1ecb58fa4497cce8e31bdaac277aed8015bbb9b9,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1fe263243300895331b3243428988dd843d4848f24b8c3ad2fa39036b85bb858,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x303297a2b9dda1f768470c971cd7d07b93fea0a7c3fe70ceb309c5ff488da096,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x29761aa4addb134d35d6c6f0dc6bea91a8addbbdd5a845d42a305aa793633f3a,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x5e1da55ef63ef460a4b5375b975ef2ff299a84a11cbb61d2e6b84fbf12a66a6,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0xc2ab09ec23c3fa4df1de8d8fd182dbb6c71cdb4ba704c5ed6bd02716e446679,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x2a46732097cec88d25191be82518307073861163ddf1a869246c2e71452dca30,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x16c5327400e4bb4fd1b67ac6a377e54b7a8706922cfe8b0e9205771e651ef2e4,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x392bc06829845aced21dda600fc4ea1e521fc2572dbc23853305258b7a13ba2,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x271371f9e82e7a1213b1200807ceeb379dff052223f458d471368b6f021a7f50,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x227a062e79bf3fc19bc3403a1ab1324f1c5dde97d32cae35bdc83bc8a1afe387,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1f5c1fff356d9f3dae47b00347a450e6021fb1eac9741278f9de13b71bd18e96,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x27121a9427d038a5c1569272ff07e4b4e4f83ff904dfe66f7f9117d121b8dccb,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x2830feb75187e962f73a7ce3454b77fb9aaf30ace72b967dc0af5c659d270427,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x9b86592c190d7d66b339a512d8be84dc81050e59f2a84021da2b95f317aae49,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x9434ff9ef157c661b12896b9c8ed24b3ecec33bee577283e542f2549de08a85,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x1422b6ec113634cc0cfe7ae92a0322926f5f80cfd51ace7a5cea760443d4f5d4,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1eb0cb83f40633bae6e92a31e1f93e7fe67a5c162aeb958776590fc42100abb1,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1e085c8cc60d318411858c26e68d8c9fb219504114766c9425cd0e90a51213e2,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2e7e88e986f81b1369fb410794f062f35bf2e174e80c03d5b8fb05507a32b0ae,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x990db0888f9f42b78f66a7066ce7909430cb6eeb2491732ba9b016415c17931,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x600aa5bdb997208ea1147d7c68cb25a36435d89186cf038464e823e090cb416,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x15fbce0379d1288e5ee663fbdc73a6948c283ab63b0669def9c8ae73a58fd013,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x269b2c2e5b355b816fb72d220fa508ac1fa2868aac2cccd72dbd60e08ec769e9,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0xe27b645febbd952392808f797bae1bf9fd11ca8e059d466174a03e7375c4c03,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x1c9f90bb5b64b0e2efecb5265ed2df9356f3b84eefa63bc2efac3fa2dd864162,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x2e258bb2029f00833075bed31af038db7d7b72c46a4c65f68d85905a177bde5c,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x824a9945071ae369b48378e9a43f43cb77bc88d1bec4bcfe9a9db7a442e3de5,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x9371025549f7dfc8dfc2bea5b64d934eed62015bee532d975a432345913d509,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x5d20472a9e23f59f7f106ab0411d2732e6afc5eb234ce3cf448d246e671edb7,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x53072023d7d2128a08a01cf66a0ac57db9f9b2baec460558662dbcf673d2bf8,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x18b09cc4ab7424e69daeb3dad5ff45ff26e3e360d41871223627cda5a81bb7dd,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0xdb4d3de74a9e21fb6495ac55ca54df84148d8d904742819afdadc60788ceec8,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x296f218296eafb40a31e644fe2d7046289c7ffd54934680dfe00be6d7ba513a6,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x21cf54260565956558de3ce7ab3e98a55ef38a9278309addb2056da4d1a93ee7,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x7addeda29cc0a4e51586f4ef5f2abdbc2a973710ebc04ad2261002a211af99c,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x2aa024b7e9f7c6d38e70e995851f7d24cf663ea24339e2e93fda1ad87094ffe7,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x2fdf4e5386682b70cb0c97984f8faba5a117a1d6686a4f27d8206d23b5e9baa1,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x29428eb9630a89c9c95c64dbac6bdff2ecfff6d8c505c09d019f4df7e8e58114,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2b8aa5cb56923bc32b865f2f178d285e71cd9f93289e4d85cf31235455491b64,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x29bda63ba967becca89898b46fd6ee4c3ebe95f86dc2b1b77af6acb243ff8720,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x25d2ea8dd1a75bd35b7f0ed5577f26b7dfd62b399cb5d63d3b46c547e05467ff,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1946c5a3abaa75544ae2596c14da4ccd63e918115d639afc975d379de0f95493,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0xe599f99097086c134289ceee8fc4bbbfba75ef42e25c9f0c6129012382c0f93,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x16d2f9d043cf918715234c7942db914891c16052c5013785e4456912086f996,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x13b0a5606183e85c0bb41afe2c7aee3af6f516c48723fb305d1253f7eed28a2c,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x27edc2f3789ef10ddb5ec6b891700a2a0aa1dd7377524e5df6a39d4d85d9be94,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1f7158235620112274cd083b4b499a7df044c6fddb9dd6ad68f1055757b2de59,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x13adefdde99a49e9c25123b443e06d56c9ebfce040ab86e8e0ada00d96be8be0,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x28c41d23719d1b36c8ece73242e99ddc3779b7cfa0cfa2c377f7d3cf8d779eb5,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x214bc6110bc38fd416fd190e5bfef3c2c5819bc109cecc2f5dfb429546d69482,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x114ff5d4e8cb82b0e55c31389119414e999196b596f7bc7da073330d212bf864,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x187b8e4ab2509e5a2ba34b0df63effd51e91c1067df5d347def17f66ef71a35c,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(
                x_den,
                mulmod(
                    0x20fd217b58f61d63d2166ce84461783d153a8bd50df402fbbc734480be2660af,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0xb45795f44d0011ce375284b58f9a5afbca608b2ab8730690efe063900ea12d0,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0xd44baf4c27986f38b5da7ab1e2e08e16a1e4516219750538324b6ef77e68874,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x2d3ee078bc1a9a94ade624f53e2c29162091a15c4d8f818de0d19b3e0fc68e11,
                    x_i,
                    N
                ),
                N
            )
            x_den := addmod(x_den, x_i, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x1f4521ac40562bd4288cab27acf0bb2d218a38547d1bcc21d988e8af252b2957,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x10075d6ca97fd780e45455ed38448f24e624581e07b0cbf2fbb00dd3da202824,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            x_num := addmod(
                x_num,
                mulmod(
                    0x155a5f113ab2d6145236199cf0ffee82df33e23ed67f6e0ab196319b5389e3,
                    x_i,
                    N
                ),
                N
            )
            y_num := addmod(
                y_num,
                mulmod(
                    0x1dbf48ac58561a63bb4090724570c67dfa795e26186754371a7e22dee2694a2,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x166eeb1df3891247d90b68201873ca954aad7f6cb3a4f52fbbb86ef9ed5c1ede,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x196937deb02818a936a48fe5add2ce3ce0bffbbe79ca150ab1d379032e101441,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1bfe16a7544597542056b8d05915e878c0e9578c1f772df5e5bff84b858b8bf7,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0xd57e8d1c71177b98c77e5d3258ee9634c95ebae53cb6ae70a2ae3d386d69a31,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2e6ed2c4b64d0f67e8e068966349ef62ad94b01ff4482c96852f398ba8f7b966,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x5b3efe5e151c00d2d22e7db0715fb72b4895453882cf7209925c32e4f28c99,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0xd307bce1fa0c2f6ab6fce9a8eecb4e2916bfa5c8cef0dbdbbbb9c11a15ba1d4,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x1d380d27775405b517c8efc5f2dd7c9b9b5ea1114da5c6a17579e8e9c167bc08,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x200d68ea21fa96861734a306ed63eeb4f90fb53ec96f214a3f8666d285b43ba1,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x1806f37b75880c0a6389731d700ba3cdae3f780fafd99433218587ce228d5c99,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1971330203ab214333667c6d6f3792369cb80c2e15481367e537476dcdda31,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0xc19759d9b6ae5368bd667a0e19a18974459e84244b2d34c2c824d06ea95dd23,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2384d43228d8f93396b569d08366e70a06a522e7347d81cadd4cd599f4b269a9,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x560b1c3b695af3b1ccbc67906904cc9477f48840c4598cb2325ec0db05f041d,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0xc99caea2136d30481b2ce433f3669f0bdd16e4b149fdfeff8cfeec4303a9585,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x2f592e2dd2ab3b7644c44cfc376803898e29de27b105631ec11ed9f6cbe65f7,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x19a2904b5c6ef16aee1af579f6c452e2a99265a9150ed44efe6c52c2d7ad7231,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0xc4ccb576d4e594a7b21980a23afdd2a7b3979cce02173f674c5f11cc46ae18b,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0xdc1da234e40cc147dd94f1e41ceae31c1c6f664ff94c9667dd603c8cf9c62d3,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x1621a9b158779ffdbefd98c8c23eb1642831c6140e72908007488ab6c2315009,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2e9020ef1af8ece6911410e48d6c8a5d17a5ef3f95e0079c2ee102c48d0dd03c,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x11fc26ce31a93ee4e07f102f63f38df46b1f383d9780e819613cb43ccb261e1a,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2358c46144014af0e9bcce0f13aae3782d3af8aababb87632592467810c7b6b3,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x178b82958c3f88deb2bb0404398397ff57f1f243f3c8d3f8655845a2fab0555b,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x10da0eb8256790bdc23046b49681a98a069f8752d7dd325d307a253622266e73,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0xceb8e0b42acf4e05ae29aab76c2968597025b82ce756e531a3680689be0e6e2,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x57274f3683f76cd7f2c56cd17ea53aab8e7e94412a080c7236d2ec670be64a,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x28e8add194b542ba1c2e9894864c672eb07b51a46bc75129e2fc187c7ade1e9e,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x6871670efe93f598743793ae822dae5d11832f9c0a755a0549adf562f826ae3,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x3402604b4b133be349c91d8115aab01e7e70318989f6ccfca06b0ccb6a4042d,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x82f33a1370a46559596bed63212a542e244f3098d55957f517de2a71a8e127e,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x2722a965620a38a089aee3b9a47c331244741611f35ad75a7b4c3f00b0d38b1c,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x10da8134b7c5b67761e83ab89fed284db77fb4cdac964ff93f4c3de298572e16,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x29c6302b2be7b3310abac2e37e2cb44027c81f9581dcb60bdf053d9125de11ef,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x22e6e91aad2027783e0105414f4b3d6278bb337442dd3af2f585973cfc9cc6a9,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x18acfff3ef41ed22c314ea6e7bc06f46ad971cb41db99a9dc43aba2e302fd75,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x24aa1251894f6fb745fd84a9ecf170b2564cf999ebf2b4e78034378ee2f2ce05,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x11f13b6c6a55a62f6f087cdd217697c30ed590e3381588df152add3b024fd29e,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x22b8fb2bb1c805de7f1e2291a70971fa46e8738e7cfdf3deafff355667858f55,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x2dec5da3fc0997e6abe182a97454dde4442f529f310516f476d8d6d8312b5b4a,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x113b8bc58f226da2eae100970b65eb9ef66e22d053b805b7f21bde987b297268,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x10ca666ea5545c131ed6acab5115c4772c12f4cc130556f86620cfaac2551bb8,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2e6df51e47a09f60b92c24ea53292397bec8b2ecbb54548933ef1159b23f35bb,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x30b580f18eda3575f906145aee24732f475174db0d52f162fdb592231f29dc5,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x1e8b5eea03d62a766dbfe3126657c25d041fbbe30352cc352af2c6df022a2fee,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x1b4227ca1b6ee41aa4432a41a7f82cc5ed18cfc7fbaa5b4423954ccf4d4b10c8,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0xee7bf087d07a0733aa37fb580880a55c7900b0dd14c9c823fab1579f6d6826f,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0xa97bb8686a74d0b441588d51ffd5e2ec273faee25033269f1c1a4d4d9b1d35,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x157cb44b51fb2714d6d7e188bbe33d86ef098429be5e1d6dabad0377fa30a634,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x24eea3c509c2ce5b910b1f37f450c3ed891b166785b742e345e8a539da5618da,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x2693bbe61812cb3cdec890036060f300f3e5f8bc73daaa1f925312a516aee340,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x1c924c263c6c95639fdd26135defcf9b9c8b1b1731ec46b522cdd0e37d7d3e52,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x23522cd01978062f007fa60a17821f4b51fcf94ef2cabfa05a47321c8833fe9,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x8f1350ce39e5ca25caccf6b36f9bf00722c10ec3354e2d2c14ba78541a1d9f4,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(
                y_den,
                mulmod(
                    0x19498aff94d85c00def9808125d1882cd4171c76e0b51f32fc9ca0b5b0fb1263,
                    x_i,
                    N
                ),
                N
            )
            x_i := mulmod(x_i, xPrime, N)
            y_num := addmod(
                y_num,
                mulmod(
                    0x1067c9af0fe2169c9e47fe83b968565c2ce5e3eb5ce072361490365f5d5fda5d,
                    x_i,
                    N
                ),
                N
            )
            y_den := addmod(y_den, x_i, N)
            x_i := mulmod(x_i, xPrime, N)
        }
        uint256 x_den_inv = inverse(x_den);
        uint256 y_den_inv = inverse(y_den);
        assembly {
            x := mulmod(x_num, x_den_inv, N)
            y := mulmod(mulmod(yPrime, y_num, N), y_den_inv, N)
        }
    }
}
