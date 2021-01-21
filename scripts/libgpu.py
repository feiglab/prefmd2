import os


def update_gpu(properties, options):
    # cuda_devices = os.getenv("CUDA_VISIBLE_DEVICES")
    if options['gpu_id'] is not None:
        if options['openmm']['platform'] in ("CUDA", "OpenCL"):
            # print("- Assigning a task to GPU: %s" % (options['gpu_id']))
            properties['DeviceIndex'] = str(options['gpu_id'])
    return properties
