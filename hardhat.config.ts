import type { HardhatUserConfig } from 'hardhat/config'
import '@nomicfoundation/hardhat-toolbox'
import 'hardhat-contract-sizer'
import 'hardhat-gas-reporter'
import * as dotenv from 'dotenv'

dotenv.config()

const config: HardhatUserConfig = {
    solidity: {
        version: '0.8.23',
        settings: {
            viaIR: false,
            optimizer: {
                enabled: true,
                runs: 100,
            },
        },
    },
    networks: {
        hardhat: {
            chainId: 1,
            blockGasLimit: 10_000_000,
            accounts: {
                count: 10,
            },
            allowUnlimitedContractSize: true,
        },
    },
    gasReporter: {
        enabled: true,
        currency: 'USD',
        gasPrice: 1,
    },
    etherscan: {
        apiKey: {
            mainnet: process.env.ETHERSCAN_API_KEY as string,
        },
    },
    contractSizer: {
        alphaSort: true,
        disambiguatePaths: false,
        runOnCompile: false,
        strict: true,
    },
}

export default config
