## Docker: build image and run/execute container

`Dockerfile` includes the instructions to build the docker image of an Ubuntu machine
with Intel-QS already installed. The image can be 'run' to create a container.
The container can be 'executed' to login into the machine.
For this, three steps need to be followed to prepare containerirzed intelqs simulator.

*  Create container
*  Configure SSH tunneling
*  Launch Jupyter notebook

#### Create container

Important: If user is not `root`, add `sudo` before each bash and docker command.  
In the code snippet below, the first command will build the docker image for Intel-QS
The default namne is qhipster, but one can choose other names too.
In the latter case, replace the chosen name in all other commands as well. 
The second command will create the container while mapping 8080 port to the localhost. 
The third command will enable conda env, execute cmake and make command, and launch Jupyter
notebook in the terminal. Copy and save the printed token from the terminal.  
`Example of a given token : http://127.0.0.1:8080/?token=6ee42173ee71353c1f1b33f8feb33132aed15f2a07960bc8`.
Last command is optional, this can be executed to loging to the container.

```bash
  docker build -t qhipster .
  docker run -d -t -p 8080:8080 qhipster
  docker exec -i $(docker ps|grep -i qhipster|cut -d ' ' -f1) /bin/bash -c ". ~/.bashrc && . /opt/intel/mkl/bin/mklvars.sh intel64 ilp64 && . /opt/intel/bin/compilervars.sh -arch intel64 -platform linux && mkdir build && cd build && CXX=g++ cmake -DIqsMPI=OFF -DIqsUtest=ON -DIqsPython=ON .. && make && cd .. && jupyter notebook --ip 0.0.0.0 --port 8080 --no-browser --allow-root"
  docker exec -it $(docker ps|grep -i qhipster|cut -d ' ' -f1) /bin/bash
```

If Docker is used on a Windows host machine, the last line should be substituted by:
`winpty docker exec -it <container_id> //bin/bash`.


#### Configure SSH tunneling

In local laptop, open Ubuntu emulator. This allows using the SSH protocol for port forwarding.
For example, if you use MobaxTerm tool, launch a session and type following command in the mobxterm shell:

```bash
  ssh -L 8080:localhost:8080 user@domain.com
```

#### Launch Jupyter notebook

Now, paste the copied token (which we copied from 'Create Container' section) in your preffered web browser 
(most importantly, please clear your browser cache before pasting the token). Once you are seeing all 
folders and files in browser, follow below section to begin.

