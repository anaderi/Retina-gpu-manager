This repository was originaly forked from the CERN repository.
https://gitlab.cern.ch/dcampora/lhcb-gpu-manager

GPU Manager
===========

The purpose of this project is in integration of accelerator hardware for
massively parallel computation with Gaudi. The GPU Manager is a process that
runs outside of Gaudi and manages execution of dedicated algorithms on the
hardware. It is capable of receiving  multiple concurrent data streams from
Gaudi and other sources.


What's inside
-------------

Brunel v45r0 is shipped with this release, alongside a modified version of the
Fast VELO Pixel tracking algorithm running on the GPU Manager. There is  a
number of supporting packages and several sample projects demonstrating how to
extend GPU Manager yourself.


Build instructions
------------------

* To compile everything:

	. init.sh
	cd Brunel_v45r0
	make

* To run the server daemon:

	. run_server.sh

* To run Brunel with the config file:

	. run_brunel_client.sh


Package Descriptions
--------------------

	Brunel_v46r0
	├── CUDA - cuda_library document type for using CUDA within CMT
	├── GpuManager
	│   ├── GpuAlgorithm     - Gaudi interface to SampleGpuHandler
	│   ├── GpuDriver        - plays recorded data back to the GPU Manager
	│   ├── GpuHandler       - base class for GPU Manager handlers
	│   ├── GpuIpc           - inter-process communication library
	│   ├── GpuServer        - main GPU Manager executable
	│   ├── GpuService       - Gaudi service communicating with GPU Manager
	│   └── SampleGpuHandler - very simple GPU Manager handler
	└── Pr
	    ├── PrPixel              - Gaudi algorithm for PrPixel
	    ├── PrPixelCuda          - GPU Manager handler for PrPixel
	    ├── PrPixelSerialization - shared serialization code
	    └── PrPixelTest          - unit tests
