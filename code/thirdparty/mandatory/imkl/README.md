The Intel Math Kernel Library is a component of many (highly-optimized) libraries by Intel.
That library suite has seen a number of changes over the years that resulted in the integration
of the IMKL library into what is called the Intel oneAPI library package. This library should be
installed via the official documentation of Intel 

```
sudo wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \ | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
sudo apt install intel-basekit
``` 
