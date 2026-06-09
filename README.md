# PNM-LIB-Cathode

Pore Network Model reproduction of:
> Khan ZA, Elkamel A, Gostick JT. "Pore Network Modelling of Galvanostatic Discharge
> Behaviour of Lithium-Ion Battery Cathodes." J. Electrochem. Soc. 168(7), 070534 (2021).

## Setup

```bash
pip install -r requirements.txt
```

## Run Tests

```bash
pytest tests/ -v
```

## Run Discharge Simulation

```bash
python scripts/run_discharge.py
```

## Project Structure

See PROMPT.md for the full design document and implementation plan.
