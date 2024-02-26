# RFC9380 Map-to-Point Implementations in Solidity

## Shallue-van de Woestijne (SVDW)

The original Shallue-van de Woestijne mapping is implemented as per [6.6.1](https://datatracker.ietf.org/doc/html/rfc9380#section-6.6.1).

Mean gas cost of SVDW is ~22.6k.

## Shallue-van de Woestijne-Ulas (SSWU)

Simplified SWU mapping is implemented as per [6.6.3](https://datatracker.ietf.org/doc/html/rfc9380#section-6.6.3) given $AB == 0$ for BN254.
We use the 59-isogeny curve $E'$ defined by:

$y'^2 = x'^3 + 9087994317191712533568698403530528306233527979934880849865820425505218365052x' + 3059101143800926337153883959975852125336293569895750485959800095292563537400$

Without the required `isoMap`, `mapToPoint` only consumes ~16k gas. Unfortunately, the 59-isogeny is the lowest-degree curve isogenous to BN254, and so `isoMap` dominates gas usage while also taking up a lot of contract bytecode.

Mean gas cost of SSWU is ~54.2k.

## ModExp

The addition chains in `ModExp.sol` are borrowed from https://github.com/thehubbleproject/hubble-contracts/blob/master/contracts/libs/ModExp.sol.

## References

See https://github.com/kevincharm/draft-irtf-cfrg-hash-to-curve/pull/1 for reference implementation Sage scripts and to generate precomputed constants.
