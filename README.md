# RFC9380 Map-to-Point Implementations in Solidity

## SSWU

Simplified SWU mapping is implemented as per [6.6.3](https://datatracker.ietf.org/doc/html/rfc9380#section-6.6.3) given $AB == 0$ for BN254.
We use the 59-isogeny curve $E'$ defined by:

$y'^2 = x'^3 + 9087994317191712533568698403530528306233527979934880849865820425505218365052x' + 3059101143800926337153883959975852125336293569895750485959800095292563537400$

Without the required `isoMap`, `mapToPoint` only consumes ~16k gas. Unfortunately, the 59-isogeny is the lowest-degree curve isogeneous to BN254, and so `isoMap` dominates gas usage while also taking up a lot of contract bytecode.
