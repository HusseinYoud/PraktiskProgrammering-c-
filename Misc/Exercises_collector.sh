#!/usr/bin/env bash


SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

cd "$SCRIPT_DIR/.." || exit 1

mkdir -v -p Exercises

for week in Week*/; do
    week_name=$(basename "$week")
    if [ -d "${week}exercises" ]; then
        for exercise in "${week}exercises/"*/; do
            [ -e "$exercise" ] || continue

            ex_name=$(basename "$exercise") #Saves the name of the exercise

           cp -r "$exercise" "Exercises/${ex_name}"
        done
    fi
done

echo "Copied all the exercises into the Exercises folder."