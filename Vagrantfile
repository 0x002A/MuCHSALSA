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

  # Disable automatic box update checking. If you disable this, then
  # boxes will only be checked for updates when the user runs
  # `vagrant box outdated`. This is not recommended.
  # config.vm.box_check_update = false

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine. In the example below,
  # accessing "localhost:8080" will access port 80 on the guest machine.
  # NOTE: This will enable public access to the opened port
  # config.vm.network "forwarded_port", guest: 80, host: 8080

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine and only allow access
  # via 127.0.0.1 to disable public access
  # config.vm.network "forwarded_port", guest: 80, host: 8080, host_ip: "127.0.0.1"

  # Create a private network, which allows host-only access to the machine
  # using a specific IP.
  # config.vm.network "private_network", ip: "192.168.33.10"

  # Create a public network, which generally matched to bridged network.
  # Bridged networks make the machine appear as another physical device on
  # your network.
  # config.vm.network "public_network"

  # Share an additional folder to the guest VM. The first argument is
  # the path on the host to the actual folder. The second argument is
  # the path on the guest to mount the folder. And the optional third
  # argument is a set of non-required options.
  # config.vm.synced_folder "../data", "/vagrant_data"

  # Set default provider.
  # config.vm.provider "virtualbox"
  # Provider-specific configuration so you can fine-tune various
  # backing providers for Vagrant. These expose provider-specific options.
  # Example for VirtualBox:
  config.vm.provider "virtualbox" do |vb|
    # # Display the VirtualBox GUI when booting the machine
    # vb.gui = true

    # Customize the amount of ressources on the VM:
    vb.memory = "4400"
    vb.cpus   = "4"
  end
  #
  # View the documentation for the provider you are using for more
  # information on available options.

  # Enable provisioning with a shell script. Additional provisioners such as
  # Ansible, Chef, Docker, Puppet and Salt are also available. Please see the
  # documentation for more information about their specific syntax and use.
  # possible deps: cmake-extras build-essential
  config.vm.provision "shell", inline: <<-SHELL
    apt-get update
    apt-get install -y clang-12 clang-format-12 clang-tidy-12 clang-tools-12 \
                       libc++-12-dev libc++abi-12-dev cmake make ninja-build \
                       doxygen graphviz \
                       jellyfish tree

    # Mambaforge (conda + mamba)setup
    echo 'export MAMBA_NO_BANNER=1' >> /etc/bash.bashrc  # kill annoying art
    echo "I am $(whoami)"
    # su - vagrant    # become user vagrant
    sudo -u vagrant bash
    echo "Now, I am $(whoami)"
    cd              # download into home directory
    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
    bash ./Mambaforge-Linux-x86_64.sh -b
    ./mambaforge/bin/conda init bash
    bash            # start new shell to load base env
    mamba update -y conda mamba

    # Install dependencies
    mamba install -y abyss python minimap2 bbmap numpy networkx biopython

    # I want this d*%&n alias ...
    echo "alias l='less -SiR'" >> ~/.bashrc
  SHELL
end
