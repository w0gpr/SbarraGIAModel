# The Sbarra model
This project is a based on the originally available code by Chris Sbarra, University at Buffalo Geology Master's Student under Jason Briner. It is a model for isostatic adjustment based on varying extents of the Greenland Ice Sheet. It is a simple spatial 1D model through time.

I am adapting this model to by used in a more dynamic way to eventually be incorporated in GHUB at the University at Buffalo.

## The Goals:
[] Allow for user defined locations (2 points).
[] The model will generate a profile along that line of present ice sheet height and bedrock surface elevation based on BedMachine V3.
[] The model will then generate the depression due to various extents and thicknesses of ice based on the bed (more details, see Chris's masters)

### The model will allow for:
[] interactive point location picking,
[] input of vector files such as continental shelf outline
[] input of ice divide
[] input of moraines and ages
[] input of database of marine limits with ages and errors 
[] weighted distance of the marine limits from the generated line

### Outputs
[] sea level curves output at selected points
[] retreat rates(?)

Development of the adaptation is currently early stages. I helped Chris a little during his development of the project and need to continue diving into the code to understand it's mechanics better. 
