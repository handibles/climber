
#0. tricks
# look for something available in repositories, good candidate for alias
apt-cache search dropbox




#1.installing R 
sudo apt-get install r-base-core

	#1.1 add self user 'user' to group 'staff' which can write to /usr/local/lib/R
	sudo adduser user staff

  #1.2 install R from source, in local dir:
  
  '[The easiest way to do this is to install R from source](https://unix.stackexchange.com/questions/149451/install-r-in-my-own-directory)':
  
  wget http://cran.rstudio.com/src/base/R-3/R-3.5.1.tar.gz
  tar xvf R-3.5.1.tar.gz
  cd R-3.5.1
  ./configure --prefix=$HOME/R
  make && make install

#2. installing RStudio 
wget https://download1.rstudio.org/rstudio-xenial-1.1.463-amd64.deb
sudo apt-get install rstudio-xenial-1.1.463-amd64.deb
rm rstudio-xenial-1.1.463-amd64.deb


#3. browser of choice - vivaldi
# get key, add key
wget http://repo.vivaldi.com/stable/linux_signing_key.pub ; sudo apt-key add linux_signing_key.pub

# add PPA to sources file - this fails:  "sudo cat '...' >> list"  ::  '>>' not as root, pipe instead
echo 'deb http://repo.vivaldi.com/stable/deb/ stable main' | sudo tee -a /etc/apt/sources.list
	

## see call down the moon
# R Stuff

