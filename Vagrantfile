# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure("2") do |config|
  # The most common configuration options are documented and commented below.
  # For a complete reference, please see the online documentation at
  # https://docs.vagrantup.com.

  # Every Vagrant development environment requires a box. You can search for
  # boxes at https://vagrantcloud.com/search.
  config.vm.box = "ubuntu/focal64"

  # Share an additional folder to the guest VM. The first argument is
  # the path on the host to the actual folder. The second argument is
  # the path on the guest to mount the folder. And the optional third
  # argument is a set of non-required options.
  # config.vm.synced_folder "../data", "/vagrant_data"

  # Set default provider.
  # config.vm.provider "virtualbox"
  # Provider-specific configuration.
  config.vm.provider "virtualbox" do |vb|
    # # Display the VirtualBox GUI when booting the machine
    # vb.gui = true

    # Customize the amount of ressources on the VM:
    vb.memory = "4400"
    vb.cpus   = "4"
  end

  # Provisioning
  config.vm.provision "root", type: "shell", privileged: true, inline: <<-SHELL
    echo "I am $(whoami)"
    echo "PWD: $(pwd)"
    apt-get update
    apt-get install -y clang-12 clang-format-12 clang-tidy-12 clang-tools-12 \
                       libc++-12-dev libc++abi-12-dev cmake make ninja-build \
                       doxygen graphviz tree

    # Global mamba config
    echo 'export MAMBA_NO_BANNER=1' >> /etc/bash.bashrc  # kill annoying art
  SHELL

  # Provisioners are executed in the order specified. The "after" option is
  # still experimental, so we don't use it here.
  config.vm.provision "user", type: "shell", privileged: false, inline: <<-SHELL
    ##### Mambaforge (conda + mamba) setup
    echo '### Installing mambaforge'
    # Download installer
    wget --no-verbose --no-clobber \
        https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
    rm -rf ./mambaforge     # clean up first (for re-provisioning)
    bash ./Mambaforge-Linux-x86_64.sh -b        # install mambaforge
    eval "$(./mambaforge/bin/conda 'shell.bash' 'hook')"  # init conda
    conda init bash                             # init conda for user

    # Set up Bioconda channels
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    # Install dependencies
    echo '### Installing dependencies using mamba'
    conda activate base
    # mamba update -y conda mamba           #  this breaks mamba o.O
    mamba install -y python numpy networkx abyss bbmap biopython \
                     kmer-jellyfish minimap2

    # I want this d*%&n alias ...
    {
        echo "alias l='less -SiR'"
        echo 'cd /vagrant                 # working dir is mounted here'
    } >> ~/.bashrc
  SHELL
end
