import { ethers } from 'hardhat'
import { MapToPointTest, MapToPointTest__factory } from '../typechain-types'
import { SignerWithAddress } from '@nomicfoundation/hardhat-ethers/signers'
import { expect } from 'chai'

describe('Map-to-point', () => {
    let deployer: SignerWithAddress
    let mtp: MapToPointTest
    beforeEach(async () => {
        ;[deployer] = await ethers.getSigners()
        mtp = await new MapToPointTest__factory(deployer).deploy()
    })

    it('correctly implements SvdW', async () => {
        const u = 7105195380181880595384217009108718366423089053558315283835256316808390512725n
        const [p0Impl, p0Gas] = await mtp.test__svdw_mapToPoint(u)

        console.log(`svdw_mapToPoint(${u}) = ${p0Impl}`)
        console.log(`p0Gas: ${p0Gas}`)
        expect(p0Impl[0]).to.eq(
            15712026073284912390314437469998384224444098668487062629391055065992760594476n,
        )
        expect(p0Impl[1]).to.eq(
            12286200326952730997678485294504458874299852441720220164574895986935631271221n,
        )
    })

    it('correctly implements SSWU', async () => {
        const u = 7105195380181880595384217009108718366423089053558315283835256316808390512725n
        const [p0Impl, p0Gas] = await mtp.test__sswu_mapToPoint(u)

        console.log(`sswu_mapToPoint(${u}) = ${p0Impl}`)
        console.log(`p0Gas: ${p0Gas}`)
        expect(p0Impl[0]).to.eq(
            7433244435151743403934667274157583038597013229141355912918907345679928483392n,
        )
        expect(p0Impl[1]).to.eq(
            3341345691842296612745507125415299735564087771630588448932624272206506288268n,
        )
    })
})
