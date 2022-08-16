import typer
from source_model import SourceModel, SourceModelAnalytical
from gui import SourceModelGui
from tools import model_several_sources, set_logger, get_analytical_image, check_if_test_source, get_image_from_path
from constants import help_analytical, help_real, help_directory


def main(terminal_anl: bool = typer.Option(False, '-a', help=help_analytical),
         source_path: str = typer.Option(None, '-r', clamp=True, help=help_real),
         dir_path: str = typer.Option(None, '-d', clamp=True, help=help_directory)
         ):
    if terminal_anl:
        set_logger(log_path='QuasarModelLog', console=True)
        image, gauss_vec = get_analytical_image()
        source = SourceModelAnalytical()
        source.gauss_vec = gauss_vec
        org, mdl, anl, _ = source.process(image)
        source.plot_all(org, mdl, anl)

    elif source_path:
        set_logger(log_path='QuasarModelLog', console=True)
        source_path = check_if_test_source(source_path)
        image = get_image_from_path(source_path)
        source = SourceModel()
        org, mdl, anl, _ = source.process(image)
        source.plot_all(org, mdl, anl)

    elif dir_path:
        set_logger(log_path='QuasarModelLog', console=True)
        model_several_sources(dir_path)

    else:
        SourceModelGui().run_gui()


if __name__ == "__main__":
    typer.run(main)

