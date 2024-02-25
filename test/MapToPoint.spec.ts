import { ethers } from 'hardhat'
import { MapToPointTest, MapToPointTest__factory } from '../typechain-types'
import { SignerWithAddress } from '@nomicfoundation/hardhat-ethers/signers'
import { expect } from 'chai'
import { hexlify } from 'ethers'
import SVDW_TEST_VECTORS from './vectors/svdw'

describe('Map-to-point', () => {
    let deployer: SignerWithAddress
    let mtp: MapToPointTest
    beforeEach(async () => {
        ;[deployer] = await ethers.getSigners()
        mtp = await new MapToPointTest__factory(deployer).deploy()
    })

    it('correctly implements SvdW', async () => {
        for (const { u, p } of SVDW_TEST_VECTORS) {
            const [pImpl] = await mtp.test__svdw_mapToPoint(u)
            expect(pImpl).to.deep.eq(p)
        }

        // fuzz gas
        let iterations = 1000n
        let sumGasCost = 0n
        for (let i = 0n; i < iterations; i++) {
            const [, gasCost] = await mtp.test__svdw_mapToPoint(pickRandomF())
            sumGasCost += gasCost
        }
        const meanGasCost = sumGasCost / iterations
        console.log(`[svdw] mean gas cost: ${meanGasCost}`)
    })

    it('correctly implements SSWU', async () => {
        const u = 7105195380181880595384217009108718366423089053558315283835256316808390512725n
        const [p0Impl] = await mtp.test__sswu_mapToPoint(u)

        expect(p0Impl[0]).to.eq(
            7433244435151743403934667274157583038597013229141355912918907345679928483392n,
        )
        expect(p0Impl[1]).to.eq(
            3341345691842296612745507125415299735564087771630588448932624272206506288268n,
        )

        // fuzz gas
        let iterations = 1000n
        let sumGasCost = 0n
        for (let i = 0n; i < iterations; i++) {
            const [, gasCost] = await mtp.test__sswu_mapToPoint(pickRandomF())
            sumGasCost += gasCost
        }
        const meanGasCost = sumGasCost / iterations
        console.log(`[sswu] mean gas cost: ${meanGasCost}`)
    })
})

/// Pick random element from BN254 F_p, accounting for modulo bias
function pickRandomF(): bigint {
    for (;;) {
        const rand32 = crypto.getRandomValues(new Uint8Array(32)) // 256-bit
        const f = BigInt(hexlify(rand32))
        if (f < 21888242871839275222246405745257275088696311157297823662689037894645226208583n) {
            return f
        }
    }
}
