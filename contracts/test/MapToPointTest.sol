// SPDX-License-Identifier: MIT
pragma solidity ^0.8;

import {SVDW} from "../SVDW.sol";
import {SSWU} from "../SSWU.sol";

contract MapToPointTest {
    SVDW immutable svdw;
    SSWU immutable sswu;

    constructor() {
        sswu = new SSWU();
        svdw = new SVDW();
    }

    function test__svdw_mapToPoint(
        uint256 u
    ) external view returns (uint256[2] memory p, uint256 gas) {
        gas = gasleft();
        p = svdw.mapToPoint(u);
        gas = gas - gasleft();
    }

    function test__sswu_mapToPoint(
        uint256 u
    ) external view returns (uint256[2] memory p, uint256 gas) {
        gas = gasleft();
        p = sswu.mapToPoint(u);
        gas = gas - gasleft();
    }
}
