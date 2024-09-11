# Client-Aided-PPML

## Building the Project

> [!IMPORTANT]
> The project should be built on Ubuntu 16.04 ~ 22.04

The [env_setup.sh](env_setup.sh) script will handle installation of build dependencies

```bash
bash env_setup.sh
```

The [build.sh](env_setup.sh) script will compile the codes and generate executable files in `build` folder.

```bash
bash build.sh all
```
## Running the Experiments

To run the experiments, go to the `build` folder, find the executable files, and use
```bash
sudo ./server 1 1025
sudo ./server 2 1025
sudo ./client 3 1025
```
in three terminals for 2 servers and the client.
