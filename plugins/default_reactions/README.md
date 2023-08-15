# Reaction template for KIMMDY

## Installation
Should get installed together with kimmdy. If you want to install it separatly: 
* Download
* `pip install -e ./`  or for development: `pip install -r requirements.txt`

## Making your own
* Implement your reaction as a subclass of kimmdy.reaction.Reaction
* Register your Reaction class in the  **[options.entry_points]** section in the setup.cfg. The name you give here must match the entry in the config.yml for Kimmdy!





