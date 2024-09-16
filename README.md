# Client-Aided-PPML

This is the implementation for [Client-Aided Privacy-Preserving Machine Learning](https://eprint.iacr.org/2024/1196) (SCN2024).

## Building the Project

> [!IMPORTANT]
> The project should be built on Ubuntu 16.04 ~ 22.04

The [env_setup.sh](env_setup.sh) script will handle installation of build dependencies. The following command will install an earlier version of [emp-toolkit](https://github.com/emp-toolkit) from [envir.zip](envir.zip): 

```bash
bash env_setup.sh
```

The [build.sh](env_setup.sh) script will compile the codes and generate executable files in the `build` subfolder:

```bash
bash build.sh all
```

You may also use the following command to build codes from a specific folder:

```bash
bash build.sh path_to_folder
```
For example, you may use `bash build.sh semi-honest_NN` to compile codes for neural networks training in the semi-honest security.

## Running the Experiments

To run the experiments, go to the `build` subfolder, find the executable files, and use the following three commands in three terminals respectively:
```bash
sudo ./server 1 1025
sudo ./server 2 1025
sudo ./client 3 1025
```
In the commands above, please replace `server` and `client` by the names of executable files for servers and client.

## Help

For any questions on building or running the library, please contact [Xinyi Shi](mailto:xinyi_shi@brown.edu) at xinyi_shi at brown dot edu.
