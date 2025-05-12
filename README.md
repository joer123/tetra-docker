# Tetra - docker

**Works on linux**

To run the system, you must find first your RTL device, for example using `lsusb`, example:

```sh
~ $ lsusb
[...]
Bus 003 Device 013: ID 0bda:2838 Realtek Semiconductor Corp. RTL2838 DVB-T
[...]
```

Here, you must write down the bus 003 and the device 013, go to the `docker-compose.yml` file and update the "devices" section with yours.

You may want to update the credentials, so update the environment variables as desired in the compose file.

Then, you are ready to run it:

```sh
docker compose up
```

First compilation may take a while.
