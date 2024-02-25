// SPDX-License-Identifier: MIT
pragma solidity ^0.8;

import {SSWU} from "../SSWU.sol";

contract MapToPointTest {
    SSWU immutable sswu;

    constructor() {
        sswu = new SSWU();
    }

    function test__sswu_mapToPoint(
        uint256 u
    ) external view returns (uint256[2] memory p, uint256 gas) {
        gas = gasleft();
        p = sswu.mapToPoint(u);
        gas = gas - gasleft();
    }
}
