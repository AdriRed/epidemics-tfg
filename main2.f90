program main
   use net_loader
   use epidemic
   use iso_fortran_env, only: int8, int32, real64
   implicit none

   integer, parameter :: dp = real64
   integer, parameter :: ik = int32
   integer, parameter :: bk = int8

   type :: simulation_config
      character(:), allocatable :: net_filename
      character(:), allocatable :: net_name
      character(:), allocatable :: output_dir
      character(:), allocatable :: batch_file
      real(dp) :: infection_rate
      real(dp) :: recovery_rate
      real(dp) :: limit_time
      integer(bk) :: model_type
      integer :: seed
      integer :: start_node
      logical :: weighted
      logical :: save_stats
      logical :: save_events
      integer :: stats_unit, stats_steps
      logical :: has_start_node
      logical :: has_batch_file
      integer :: events_unit
   end type simulation_config


   type(epidemic_net) :: network
   type(simulation_config) :: configuration

   ! Procesar argumentos de línea de comandos
   configuration = parse_arguments()

   ! Mostrar configuración
   call print_config(configuration)

   ! Cargar la red
   network = load_network(configuration%net_filename, configuration%weighted)

   ! Ejecutar simulación
   if (configuration%has_batch_file) then
      call run_batch_simulations(network, configuration)
   else
      call run_simulation(network, configuration)
   end if

contains

   !===============================================================================
   ! Tipos de datos para configuración
   !===============================================================================


   !===============================================================================
   ! Procesamiento de argumentos
   !===============================================================================

   function parse_arguments() result(config)
      type(simulation_config) :: config
      integer :: i, n_args
      character(len=256) :: arg
      character(len=:), allocatable :: trimmed_arg  ! <-- Añadir esta variable
      logical :: infection_rate_set, recovery_rate_set, model_set, net_file_set
      logical :: help_requested

      ! Inicializar valores por defecto
      call set_default_config(config, infection_rate_set, recovery_rate_set, &
         model_set, net_file_set)
      help_requested = .false.

      ! Procesar argumentos
      n_args = command_argument_count()
      i = 1
      do while (i <= n_args)
         call get_command_argument(i, arg)
         trimmed_arg = trim(arg)

         if (trimmed_arg == '--help' .or. trimmed_arg == '-h') then
            call print_help()
            stop 0
         else if (trimmed_arg == '--infection-rate' .or. trimmed_arg == '-i') then
            call read_real_arg(i, n_args, 'infection-rate', config%infection_rate)
            infection_rate_set = .true.
         else if (trimmed_arg == '--recovery-rate' .or. trimmed_arg == '-r') then
            call read_real_arg(i, n_args, 'recovery-rate', config%recovery_rate)
            recovery_rate_set = .true.
         else if (trimmed_arg == '--start-node' .or. trimmed_arg == '-sn') then
            call read_integer_arg(i, n_args, 'start-node', config%start_node)
            config%has_start_node = .true.
         else if (trimmed_arg == '--model' .or. trimmed_arg == '-m') then
            call read_model_arg(i, n_args, config%model_type)
            model_set = .true.
         else if (trimmed_arg == '--limit-time' .or. trimmed_arg == '-lt') then
            call read_real_arg(i, n_args, 'limit-time', config%limit_time)
         else if (trimmed_arg == '--weighted' .or. trimmed_arg == '-w') then
            config%weighted = .true.
         else if (trimmed_arg == '--seed' .or. trimmed_arg == '-s') then
            call read_integer_arg(i, n_args, 'seed', config%seed)
         else if (trimmed_arg == '--events' .or. trimmed_arg == '-ev') then
            config%save_events = .true.
         else if (trimmed_arg == '--stats' .or. trimmed_arg == '-st') then
            config%save_stats = .true.
         else if (trimmed_arg == '--output-dir' .or. trimmed_arg == '-o') then  ! <-- Nueva opción
            call read_string_arg(i, n_args, 'output-dir', config%output_dir)
         else if (trimmed_arg == '--batch-file' .or. trimmed_arg == '-b') then
            call read_string_arg(i, n_args, 'batch-file', config%batch_file)
            config%has_batch_file = .true.
         else if (trimmed_arg == '--stats-steps' .or. trimmed_arg == '-ss') then
            call read_integer_arg(i, n_args, 'stats-steps', config%stats_steps)
         else
            call handle_network_file(arg, config%net_filename, net_file_set)
         end if
         i = i + 1
      end do

      ! Si no hay argumentos, mostrar ayuda
      if (n_args == 0) then
         call print_help()
         stop 0
      end if

      ! Validar argumentos obligatorios
      call validate_required_args(infection_rate_set, recovery_rate_set, &
         model_set, net_file_set, config%has_batch_file, config%stats_steps)

      ! Configurar seed si es necesario
      call setup_seed(config%seed)

      ! Obtener nombre de la red sin extensión
      config%net_name = extract_network_name(config%net_filename)

   end function parse_arguments

   subroutine print_help()
      implicit none

      write(*, '(A)') '================================================================================'
      write(*, '(A)') 'SIMULADOR DE EPIDEMIAS EN REDES'
      write(*, '(A)') 'Adrià Rojo'
      write(*, '(A)') '================================================================================'
      write(*, '(A)') ''
      write(*, '(A)') 'Uso: programa [OPCIONES] ARCHIVO_RED'
      write(*, '(A)') ''
      write(*, '(A)') 'ARGUMENTOS OBLIGATORIOS:'
      write(*, '(A)') '  ARCHIVO_RED                    Archivo de la red (formato compatible con net_loader)'
      write(*, '(A)') ''
      write(*, '(A)') 'OPCIONES:'
      write(*, '(A)') '  -h, --help                     Muestra esta ayuda'
      write(*, '(A)') ''
      write(*, '(A)') '  -i, --infection-rate VALOR     Tasa de infección (obligatorio)'
      write(*, '(A)') '  -r, --recovery-rate VALOR      Tasa de recuperación (obligatorio)'
      write(*, '(A)') '  -m, --model {SIR,SIS}          Modelo epidémico (obligatorio)'
      write(*, '(A)') ''
      write(*, '(A)') '  -lt, --limit-time VALOR        Tiempo máximo de simulación (default: 50.0)'
      write(*, '(A)') '  -s, --seed VALOR               Semilla para el generador aleatorio'
      write(*, '(A)') '  -sn, --start-node              Indica el nodo inicial para infectar'
      write(*, '(A)') '                                 (default: el nodo con degree más alto)'
      write(*, '(A)') '  -w, --weighted                 Indica que la red es ponderada'
      write(*, '(A)') '  -b, --batch-file ARCHIVO       Archivo con lista de simulaciones a ejecutar'
      write(*, '(A)') '                                 (cada línea: inf_rate rec_rate seed limit_time model [start_node])'
      write(*, '(A)') ''
      write(*, '(A)') 'ARCHIVOS DE SALIDA:'
      write(*, '(A)') '  -o, --output                   Carpeta de salida de archivos'
      write(*, '(A)') '                                 (./output)'
      write(*, '(A)') '  -st, --stats                   Guardar archivo de estadísticas'
      write(*, '(A)') '                                 (./output/stats-NOMBRE_red-I...dat)'
      write(*, '(A)') '  -ss, --stats-steps             Estadísticas cada x pasos'
      write(*, '(A)') '  -ev, --events                   Guardar archivo de eventos'
      write(*, '(A)') '                                 (./output/events-NOMBRE_red-I...dat)'
      write(*, '(A)') ''
      write(*, '(A)') 'EJEMPLOS:'
      write(*, '(A)') '  # Simulación básica SIR'
      write(*, '(A)') '  programa -i 2.0 -r 1.0 -m SIR red.txt'
      write(*, '(A)') ''
      write(*, '(A)') '  # Simulación SIS con todas las opciones'
      write(*, '(A)') '  programa --infection-rate 2.0 --recovery-rate 1.0 \'
      write(*, '(A)') '           --model SIS --limit-time 500 --weighted \'
      write(*, '(A)') '           --seed 12345 --stats --events red.txt'
      write(*, '(A)') ''
      write(*, '(A)') '  # Usar nombres cortos'
      write(*, '(A)') '  programa -i 2.0 -r 1.0 -m SIR -lt 100 -s 42 red.txt'
      write(*, '(A)') ''
      write(*, '(A)') 'ARCHIVOS DE SALIDA:'
      write(*, '(A)') '  Los archivos se guardan en ./output/ con el formato:'
      write(*, '(A)') '  stats-NOMBRE_red-I{infection_rate}-R{recovery_rate}-S{seed}.dat'
      write(*, '(A)') '  events-NOMBRE_red-I{infection_rate}-R{recovery_rate}-S{seed}.dat'
      write(*, '(A)') ''
      write(*, '(A)') '  Para redes ponderadas se añade prefijo "w":'
      write(*, '(A)') '  stats-wNOMBRE_red-I{infection_rate}-R{recovery_rate}-S{seed}.dat'
      write(*, '(A)') '================================================================================'
   end subroutine print_help

   subroutine set_default_config(config, infection_rate_set, recovery_rate_set, &
      model_set, net_file_set)
      type(simulation_config), intent(out) :: config
      logical, intent(out) :: infection_rate_set, recovery_rate_set
      logical, intent(out) :: model_set, net_file_set

      infection_rate_set = .false.
      recovery_rate_set = .false.
      model_set = .false.
      net_file_set = .false.

      config%weighted = .false.
      config%save_stats = .false.
      config%save_events = .false.
      config%has_start_node = .false.
      config%limit_time = 50.0_dp
      config%seed = 0
      config%stats_steps = 100
      config%start_node = -1
      config%stats_unit = 21
      config%events_unit = 22
      config%output_dir = './output'  ! <-- Valor por defecto
      config%batch_file = ''
      config%has_batch_file = .false.
   end subroutine set_default_config

   subroutine read_real_arg(i, n_args, arg_name, value)
      integer, intent(inout) :: i
      integer, intent(in) :: n_args
      character(*), intent(in) :: arg_name
      real(dp), intent(out) :: value
      character(len=256) :: arg
      integer :: io

      i = i + 1
      if (i > n_args) call argument_error('Falta valor para ' // arg_name)
      call get_command_argument(i, arg)
      read(arg, *, iostat=io) value
      if (io /= 0) call argument_error('Valor inválido para ' // arg_name)
   end subroutine read_real_arg

   subroutine read_string_arg(i, n_args, arg_name, value)
      integer, intent(inout) :: i
      integer, intent(in) :: n_args
      character(*), intent(in) :: arg_name
      character(:), allocatable, intent(out) :: value
      character(len=256) :: arg
      integer :: len_arg

      i = i + 1
      if (i > n_args) then
         write(*,*) 'Error: Falta valor para ', arg_name
         stop 1
      end if

      call get_command_argument(i, arg)
      len_arg = len_trim(arg)

      if (len_arg == 0) then
         write(*,*) 'Error: Valor vacío para ', arg_name
         stop 1
      end if

      allocate(character(len_arg) :: value)
      value = arg(1:len_arg)
   end subroutine read_string_arg

   subroutine read_integer_arg(i, n_args, arg_name, value)
      integer, intent(inout) :: i
      integer, intent(in) :: n_args
      character(*), intent(in) :: arg_name
      integer, intent(out) :: value
      character(len=256) :: arg
      integer :: io

      i = i + 1
      if (i > n_args) call argument_error('Falta valor para ' // arg_name)
      call get_command_argument(i, arg)
      read(arg, *, iostat=io) value
      if (io /= 0) call argument_error('Valor inválido para ' // arg_name)
   end subroutine read_integer_arg

   subroutine read_model_arg(i, n_args, model_type)
      integer, intent(inout) :: i
      integer, intent(in) :: n_args
      integer(bk), intent(out) :: model_type
      character(len=256) :: arg

      i = i + 1
      if (i > n_args) call argument_error('Falta valor para model (SIR/SIS)')
      call get_command_argument(i, arg)

      select case(trim(arg))
       case('SIR')
         model_type = SIR_MODEL
       case('SIS')
         model_type = SIS_MODEL
       case default
         write(*,*) 'Error: Modelo debe ser SIR o SIS, se obtuvo: ', trim(arg)
         stop 1
      end select
   end subroutine read_model_arg

   subroutine handle_network_file(arg, net_filename, net_file_set)
      character(*), intent(in) :: arg
      character(:), allocatable, intent(out) :: net_filename
      logical, intent(inout) :: net_file_set
      integer :: len_arg

      if (.not. net_file_set) then
         ! Obtener longitud real sin espacios al final
         len_arg = len_trim(arg)

         ! Asignar solo la longitud necesaria
         allocate(character(len_arg) :: net_filename)
         net_filename = arg(1:len_arg)

         net_file_set = .true.
      else
         write(*,*) 'Error: Múltiples archivos de red especificados: ', trim(arg)
         stop 1
      end if
   end subroutine handle_network_file

   subroutine validate_required_args(infection_rate_set, recovery_rate_set, &
      model_set, net_file_set, batch_file_present, stats_steps)
      logical, intent(in) :: infection_rate_set, recovery_rate_set
      logical, intent(in) :: model_set, net_file_set
      logical, intent(in) :: batch_file_present   ! <-- NUEVO argumento
      integer, intent(in) :: stats_steps

      if (batch_file_present) then
         ! En modo batch solo se necesita el archivo de red
         if (.not. net_file_set) then
            call argument_error('Debe especificar el archivo de red')
         end if
         return
      end if

      if (stats_steps <= 0) then
         call argument_error('--stats-steps o -ss debe ser mayor que 1')
      end if

      if (.not. infection_rate_set) then
         call argument_error('Debe especificar --infection-rate o -i')
      end if

      if (.not. recovery_rate_set) then
         call argument_error('Debe especificar --recovery-rate o -r')
      end if

      if (.not. model_set) then
         call argument_error('Debe especificar --model o -m (SIR/SIS)')
      end if

      if (.not. net_file_set) then
         call argument_error('Debe especificar el archivo de red')
      end if
   end subroutine validate_required_args

   subroutine argument_error(message)
      character(*), intent(in) :: message
      write(*,*) 'Error: ', message
      stop 1
   end subroutine argument_error

   subroutine setup_seed(seed)
      integer, intent(inout) :: seed
      integer :: seed_value

      if (seed == 0) then
         call system_clock(count=seed_value)
         seed = seed_value
      end if
   end subroutine setup_seed

   !===============================================================================
   ! Utilidades de red
   !===============================================================================

   function extract_network_name(filename) result(name)
      character(*), intent(in) :: filename
      character(:), allocatable :: name
      integer :: i, j, start_pos

      ! Encontrar el último separador de directorio
      i = index(filename, '/', back=.true.)
      j = index(filename, '\', back=.true.)
      start_pos = max(i, j) + 1
      if (start_pos == 1) start_pos = 1

      ! Extraer nombre sin extensión
      name = trim(filename(start_pos:))
      j = index(name, '.', back=.true.)
      if (j > 0) then
         name = name(:j-1)
      end if
   end function extract_network_name

   function load_network(filename, weighted) result(net)
      character(*), intent(in) :: filename
      logical, intent(in) :: weighted
      type(epidemic_net) :: net
      integer :: io, unit

      unit = 11
      open(unit=unit, file=filename, action='read', status='old', iostat=io)
      if (io /= 0) then
         write(*,*) 'Error: No se puede abrir el archivo ', filename
         stop 1
      end if

      net = initialize_net(unit, weighted=weighted)
      ! call net%hashmap%clear()
      call net%print_stats()
      close(unit=unit)
   end function load_network

   !===============================================================================
   ! Utilidades de simulación
   !===============================================================================

   subroutine print_config(config)
      type(simulation_config), intent(in) :: config

      write(*,*) '=== Configuración de la simulación ==='
      write(*,*) 'Archivo de red: ', config%net_filename
      write(*,*) 'Nombre de red: ', config%net_name
      write(*,*) 'Modelo: ', merge('SIR', 'SIS', config%model_type == SIR_MODEL)
      write(*,*) 'Infection rate: ', config%infection_rate
      write(*,*) 'Recovery rate: ', config%recovery_rate
      write(*,*) 'Limit time: ', config%limit_time
      write(*,*) 'Carpeta de output: ', config%output_dir
      if (config%has_start_node) then
         write(*,*) 'Nodo inicial: ', config%start_node
      else
         write(*,*) 'Nodo inicial: auto'
      end if
      write(*,*) 'Seed: ', config%seed
      write(*,*) 'Weighted: ', config%weighted
      write(*,*) 'Guardar stats: ', config%save_stats
      write(*,*) 'Guardar events: ', config%save_events
      write(*,*) '======================================='
   end subroutine print_config

   subroutine run_simulation(net, config)
      type(epidemic_net), intent(inout) :: net
      type(simulation_config), intent(in) :: config

      if (config%save_stats .and. config%save_events) then
         if (config%has_start_node) then
            call execute_simulation(net, config%infection_rate, config%recovery_rate, &
               int(config%seed, ik), config%limit_time, &
               config%model_type, config%output_dir, start_node=config%start_node, stats_unit=config%stats_unit, &
               events_unit=config%events_unit, net_name=config%net_name, stats_steps=config%stats_steps)
         else
            call execute_simulation(net, config%infection_rate, config%recovery_rate, &
               int(config%seed, ik), config%limit_time, &
               config%model_type, config%output_dir, stats_unit=config%stats_unit, &
               events_unit=config%events_unit, net_name=config%net_name, stats_steps=config%stats_steps)
         end if

      else if (config%save_stats) then
         if (config%has_start_node) then
            call execute_simulation(net, config%infection_rate, config%recovery_rate, &
               int(config%seed, ik), config%limit_time, &
               config%model_type, config%output_dir, start_node=config%start_node, stats_unit=config%stats_unit, &
               net_name=config%net_name, stats_steps=config%stats_steps)
         else
            call execute_simulation(net, config%infection_rate, config%recovery_rate, &
               int(config%seed, ik), config%limit_time, &
               config%model_type, config%output_dir, stats_unit=config%stats_unit, &
               net_name=config%net_name, stats_steps=config%stats_steps)
         end if
      else if (config%save_events) then
         if (config%has_start_node) then
            call execute_simulation(net, config%infection_rate, config%recovery_rate, &
               int(config%seed, ik), config%limit_time, &
               config%model_type, config%output_dir, start_node=config%start_node, events_unit=config%events_unit, &
               net_name=config%net_name, stats_steps=config%stats_steps)
         else
            call execute_simulation(net, config%infection_rate, config%recovery_rate, &
               int(config%seed, ik), config%limit_time, &
               config%model_type, config%output_dir, events_unit=config%events_unit, &
               net_name=config%net_name, stats_steps=config%stats_steps)
         end if
      else
         if (config%has_start_node) then
            call execute_simulation(net, config%infection_rate, config%recovery_rate, &
               int(config%seed, ik), config%limit_time, &
               config%model_type, config%output_dir, start_node=config%start_node, net_name=config%net_name, &
               stats_steps=config%stats_steps)
         else
            call execute_simulation(net, config%infection_rate, config%recovery_rate, &
               int(config%seed, ik), config%limit_time, &
               config%model_type, config%output_dir, net_name=config%net_name, stats_steps=config%stats_steps)
         end if
      end if
   end subroutine run_simulation

   subroutine run_batch_simulations(net, config)
      use omp_lib          ! <-- Añadir este use si no lo tienes global
      type(epidemic_net), intent(inout) :: net
      type(simulation_config), intent(in) :: config
      integer :: unit, ios, line_num, num_tasks, i
      character(len=256) :: line
      logical :: save_stats, save_events
      integer, parameter :: MAX_TASKS = 100000   ! Ajusta según necesidad
      real(dp), allocatable :: task_inf_rate(:), task_rec_rate(:), task_limit_time(:)
      integer(ik), allocatable :: task_seed(:), task_start_node(:)
      integer(bk), allocatable :: task_model_type(:)
      integer :: stats_u_base, events_u_base

      save_stats = config%save_stats
      save_events = config%save_events
      stats_u_base = 1000
      events_u_base = 2000

      ! Primero contar cuántas tareas hay
      open(newunit=unit, file=config%batch_file, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'Error abriendo archivo batch: ', config%batch_file
         stop 1
      end if

      num_tasks = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
         num_tasks = num_tasks + 1
      end do
      rewind(unit)

      if (num_tasks == 0) then
         write(*,*) 'Archivo batch vacío o sin tareas válidas.'
         close(unit)
         return
      end if

      ! Reservar memoria para las tareas
      allocate(task_inf_rate(num_tasks), task_rec_rate(num_tasks), &
         task_limit_time(num_tasks), task_seed(num_tasks), &
         task_model_type(num_tasks), task_start_node(num_tasks))

      ! Leer y almacenar tareas
      line_num = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (line(1:1) == '#' .or. len_trim(line) == 0) cycle

         line_num = line_num + 1
         ! Leer valores obligatorios
         read(line, *, iostat=ios) task_inf_rate(line_num), task_rec_rate(line_num), &
            task_seed(line_num), task_limit_time(line_num), &
            task_model_type(line_num)
         if (ios /= 0) then
            write(*,*) 'Error en línea ', line_num, ': formato incorrecto.'
            task_start_node(line_num) = -1
            cycle
         end if

         ! Intentar leer start_node opcional
         task_start_node(line_num) = -1
         read(line, *, iostat=ios) task_inf_rate(line_num), task_rec_rate(line_num), &
            task_seed(line_num), task_limit_time(line_num), &
            task_model_type(line_num), task_start_node(line_num)
         if (ios /= 0) task_start_node(line_num) = -1   ! no proporcionado
      end do
      close(unit)

      write(*,*) 'Ejecutando ', num_tasks, ' simulaciones en paralelo con OpenMP.'

      ! Bucle paralelo sobre las tareas
      !$omp parallel do private(i) schedule(dynamic)
      do i = 1, num_tasks
         if (task_start_node(i) >= 0) then
            if (save_stats .and. save_events) then
               call execute_simulation(net, task_inf_rate(i), task_rec_rate(i), &
                  task_seed(i), task_limit_time(i), task_model_type(i), &
                  config%output_dir, start_node=task_start_node(i), &
                  stats_unit=stats_u_base+i, events_unit=events_u_base+i, &
                  net_name=config%net_name, stats_steps=config%stats_steps)
            else if (save_stats) then
               call execute_simulation(net, task_inf_rate(i), task_rec_rate(i), &
                  task_seed(i), task_limit_time(i), task_model_type(i), &
                  config%output_dir, start_node=task_start_node(i), &
                  stats_unit=stats_u_base+i, net_name=config%net_name, stats_steps=config%stats_steps)
            else if (save_events) then
               call execute_simulation(net, task_inf_rate(i), task_rec_rate(i), &
                  task_seed(i), task_limit_time(i), task_model_type(i), &
                  config%output_dir, start_node=task_start_node(i), &
                  events_unit=events_u_base+i, net_name=config%net_name, stats_steps=config%stats_steps)
            else
               call execute_simulation(net, task_inf_rate(i), task_rec_rate(i), &
                  task_seed(i), task_limit_time(i), task_model_type(i), &
                  config%output_dir, start_node=task_start_node(i), &
                  net_name=config%net_name, stats_steps=config%stats_steps)
            end if
         else
            if (save_stats .and. save_events) then
               call execute_simulation(net, task_inf_rate(i), task_rec_rate(i), &
                  task_seed(i), task_limit_time(i), task_model_type(i), &
                  config%output_dir, stats_unit=stats_u_base+i, &
                  events_unit=events_u_base+i, net_name=config%net_name, stats_steps=config%stats_steps)
            else if (save_stats) then
               call execute_simulation(net, task_inf_rate(i), task_rec_rate(i), &
                  task_seed(i), task_limit_time(i), task_model_type(i), &
                  config%output_dir, stats_unit=stats_u_base+i, &
                  net_name=config%net_name, stats_steps=config%stats_steps)
            else if (save_events) then
               call execute_simulation(net, task_inf_rate(i), task_rec_rate(i), &
                  task_seed(i), task_limit_time(i), task_model_type(i), &
                  config%output_dir, events_unit=events_u_base+i, &
                  net_name=config%net_name, stats_steps=config%stats_steps)
            else
               call execute_simulation(net, task_inf_rate(i), task_rec_rate(i), &
                  task_seed(i), task_limit_time(i), task_model_type(i), &
                  config%output_dir, net_name=config%net_name, stats_steps=config%stats_steps)
            end if
         end if
      end do
      !$omp end parallel do

      deallocate(task_inf_rate, task_rec_rate, task_limit_time, task_seed, &
         task_model_type, task_start_node)

      write(*,*) 'Batch completado: ', num_tasks, ' simulaciones ejecutadas.'
   end subroutine run_batch_simulations

   !===============================================================================
   ! Función para encontrar nodo de grado máximo
   !===============================================================================

   function find_max_degree_node(net) result(node_index)
      type(epidemic_net), intent(in) :: net
      integer(ik) :: node_index
      integer(ik) :: i, max_degree

      max_degree = 0
      node_index = 1

      do i = 1, net%stats%nodes_count
         if (net%degree(i) > max_degree) then
            max_degree = net%degree(i)
            node_index = i
         end if
      end do
   end function find_max_degree_node

   !===============================================================================
   ! Generación de nombre de archivo
   !===============================================================================

   function generate_filename(net, infection_rate, recovery_rate, seed, &
      model_type, start_node, net_name) result(filename)
      type(epidemic_net), intent(in) :: net
      real(dp), intent(in) :: infection_rate, recovery_rate
      integer(ik), intent(in) :: seed, start_node
      integer(bk), intent(in) :: model_type
      character(*), intent(in), optional :: net_name
      character(:), allocatable :: filename

      character(100) :: name1, name2, name3
      ! Parte 1: parámetros numéricos
      write(name1, '(A,F10.5,A,F10.5,A,I5,A,I5.5)') 'I=', infection_rate, &
         '-R=', recovery_rate, '-S=', seed, '-SN=', start_node

      ! Parte 2: modelo
      if (model_type == SIR_MODEL) then
         write(name2, '(A, A)') 'SIR-', trim(name1)
      else
         write(name2, '(A, A)') 'SIS-', trim(name1)
      end if

      ! Parte 3: weighted
      if (net%weighted) then
         write(name3, '(A, A)') 'w', trim(name2)
      else
         write(name3, '(A)') trim(name2)
      end if

      ! Nombre completo
      if (present(net_name)) then
         filename = trim(net_name) // '-' // trim(name3)
      else
         filename = trim(name3)
      end if
   end function generate_filename

   !===============================================================================
   ! Escritura de headers
   !===============================================================================

   subroutine write_stats_common(unit, net, model_type, seed, infection_rate, &
      recovery_rate, start_node_id, start_node_degree, net_name)
      integer(ik), intent(in) :: unit
      type(epidemic_net), intent(in) :: net
      integer(bk), intent(in) :: model_type
      integer(ik), intent(in) :: seed
      real(dp), intent(in) :: infection_rate, recovery_rate
      integer(ik), intent(in) :: start_node_id, start_node_degree
      character(*), intent(in), optional :: net_name

      write(unit, *) '# Epidemic simulation - stats'
      write(unit, *) '# Weighted: ', merge('YES', 'NO ', net%weighted)
      write(unit, *) '# Model: ', merge('SIR', 'SIS', model_type == SIR_MODEL)

      if (present(net_name)) then
         write(unit, *) '# Net name: ', net_name
      end if

      write(unit, *) '# Using mt19993 random generator with seed: ', seed
      write(unit, *) '# Infection rate: ', infection_rate
      write(unit, *) '# Recovery rate: ', recovery_rate
      write(unit, *) '# Start node: ', start_node_id
      write(unit, *) '# Start node degree: ', start_node_degree
      write(unit, *) '# -----------------------------------------------------'
   end subroutine write_stats_common

   subroutine write_stats_header(unit, net, model_type, seed, infection_rate, &
      recovery_rate, start_node_id, start_node_degree, net_name)
      integer(ik), intent(in) :: unit
      type(epidemic_net), intent(in) :: net
      integer(bk), intent(in) :: model_type
      integer(ik), intent(in) :: seed
      real(dp), intent(in) :: infection_rate, recovery_rate
      integer(ik), intent(in) :: start_node_id, start_node_degree
      character(*), intent(in), optional :: net_name
      call write_stats_common(unit, net, model_type, seed, infection_rate, recovery_rate, &
         start_node_id, start_node_degree, net_name)
      write(unit, *) '# time, infected_density, recovered_density, actual_infection_rate, actual_recovery_rate'
   end subroutine write_stats_header

   subroutine write_events_header(unit, net, model_type, seed, infection_rate, &
      recovery_rate, start_node_id, start_node_degree, net_name)
      integer(ik), intent(in) :: unit
      type(epidemic_net), intent(in) :: net
      integer(bk), intent(in) :: model_type
      integer(ik), intent(in) :: seed
      real(dp), intent(in) :: infection_rate, recovery_rate
      integer(ik), intent(in) :: start_node_id, start_node_degree
      character(*), intent(in), optional :: net_name

      call write_stats_common(unit, net, model_type, seed, infection_rate, recovery_rate, &
         start_node_id, start_node_degree, net_name)
      write(unit, *) '# time, node_id, event'
   end subroutine write_events_header

   !===============================================================================
   ! Subrutina principal de simulación (adaptada)
   !===============================================================================

   subroutine execute_simulation(initialized_net, infection_rate, recovery_rate, &
      seed, limit_time, model_type, output_dir, start_node, stats_unit, events_unit, net_name, stats_steps)
      implicit none
      class(epidemic_net), intent(inout) :: initialized_net
      real(dp), intent(in) :: infection_rate, recovery_rate, limit_time
      integer(ik), intent(in) :: seed
      integer(ik), intent(in), optional :: stats_unit, events_unit, start_node, stats_steps
      integer(bk), intent(in) :: model_type
      character(:), intent(in), allocatable, optional :: net_name

      type(epidemic_simulation) :: sim
      logical :: should_write_stats, should_write_events
      integer(ik) :: start_node_index, start_node_degree, start_node_id
      character(:), allocatable :: filename
      character(*), intent(in) :: output_dir  ! <-- Nuevo parámetro

      ! Encontrar nodo de grado máximo para iniciar
      if (present(start_node)) then
         call initialized_net%hashmap%get(start_node, start_node_index)
      else
         start_node_index = find_max_degree_node(initialized_net)
      end if
      start_node_degree = initialized_net%degree(start_node_index)

      call initialized_net%rev_hashmap%get(start_node_index, start_node_id)

      ! Generar nombre de archivo
      filename = generate_filename(initialized_net, infection_rate, recovery_rate, &
         seed, model_type, start_node, net_name)

      !$omp critical(name_write)
      write(*, '(A, A)') 'Filename will be ', filename
      !$omp end critical(name_write)

      ! Inicializar simulación
      sim = initialize_simulation(initialized_net, seed, model_type)
      sim%infection_rate = infection_rate
      sim%recovery_rate = recovery_rate
      call sim%set_infected_node(start_node_index)

      ! Configurar escritura
      should_write_stats = present(stats_unit)
      should_write_events = present(events_unit)

      ! Abrir archivos y escribir headers
      call open_output_files(initialized_net, stats_unit, events_unit, filename, &
         model_type, seed, infection_rate, recovery_rate, &
         start_node_id, start_node_degree, net_name, &
         should_write_stats, should_write_events, output_dir)


      if (should_write_events) then
         write(events_unit, "(E20.10, I10, A2)") 0., start_node_id, 'I'
      end if

      ! Ejecutar simulación
      call run_simulation_loop(sim, initialized_net, stats_unit, events_unit, &
         infection_rate, recovery_rate, seed, limit_time, &
         should_write_stats, stats_steps, should_write_events)

      ! Limpiar
      call close_output_files(stats_unit, events_unit, should_write_stats, should_write_events)
      call sim%clear()

   end subroutine execute_simulation

   subroutine open_output_files(net, stats_unit, events_unit, filename, &
      model_type, seed, infection_rate, recovery_rate, &
      start_node_id, start_node_degree, net_name, &
      should_write_stats, should_write_events, output_dir)
      type(epidemic_net), intent(in) :: net
      integer(ik), intent(in), optional :: stats_unit, events_unit
      character(*), intent(in) :: filename
      integer(bk), intent(in) :: model_type
      integer(ik), intent(in) :: seed
      real(dp), intent(in) :: infection_rate, recovery_rate
      integer(ik), intent(in) :: start_node_id, start_node_degree
      character(*), intent(in), optional :: net_name
      logical, intent(in) :: should_write_stats, should_write_events
      character(*), intent(in) :: output_dir  ! <-- Nuevo parámetro

      !$omp critical(file_write)
      if (should_write_stats) then
         open(unit=stats_unit, file=output_dir // '/stats-' // filename // '.dat', action='write')
         call write_stats_header(stats_unit, net, model_type, seed, infection_rate, &
            recovery_rate, start_node_id, start_node_degree, net_name)
      end if

      if (should_write_events) then
         open(unit=events_unit, file=output_dir // '/events-' // filename // '.dat', action='write')
         call write_events_header(events_unit, net, model_type, seed, infection_rate, &
            recovery_rate, start_node_id, start_node_degree, net_name)
      end if
      !$omp end critical(file_write)

      write(*, '(A,F10.5,A,I10,A)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=start'
   end subroutine open_output_files

   subroutine run_simulation_loop(sim, net, stats_unit, events_unit, &
      infection_rate, recovery_rate, seed, limit_time, &
      should_write_stats, stats_steps, should_write_events)
      type(epidemic_simulation), intent(inout) :: sim
      type(epidemic_net), intent(inout) :: net
      integer(ik), intent(in), optional :: stats_unit, events_unit, stats_steps
      real(dp), intent(in) :: infection_rate, recovery_rate, limit_time
      integer(ik), intent(in) :: seed
      logical, intent(in) :: should_write_stats, should_write_events

      type(epidemic_step_event) :: sim_event
      type(epidemic_simulation_stats) :: sim_stats
      integer(ik) :: node_id, step_n
      step_n = 1

      do
         sim_event = sim%act()

         if (should_write_stats) then
            if (mod(step_n, stats_steps) == 0) then
               sim_stats = sim%get_stats()
               !$omp critical(file_write)
               write(stats_unit, "(E20.10, E20.10, E20.10, E20.10, E20.10)") sim%time, &
                  sim_stats%infected_density, sim_stats%recovered_density, &
                  sim_stats%rates%actual_infection_rate, sim_stats%rates%actual_recovery_rate
               !$omp end critical(file_write)
            end if
         end if

         if (should_write_events) then
            call net%rev_hashmap%get(sim_event%selected_node, node_id)
            !$omp critical(file_write)
            write(events_unit, "(E20.10, I10, A2)") sim%time, node_id, sim_event%action
            !$omp end critical(file_write)
         end if

         ! Verificar condiciones de terminación
         if (sim_event%action == 'E') then
            write(*, '(A,F10.5,A,I10,A)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=dead'
            exit
         end if

         if (sim%time > limit_time) then
            write(*, '(A,F10.5,A,I10,A)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=max'
            exit
         end if
         ! call sim%verify_consistency()
         step_n = step_n+1
      end do

      if (should_write_stats) then
         if (mod(step_n, 1000) /= 0) then
            sim_stats = sim%get_stats()
            !$omp critical(file_write)
            write(stats_unit, "(E20.10, E20.10, E20.10, E20.10, E20.10)") sim%time, &
               sim_stats%infected_density, sim_stats%recovered_density, &
               sim_stats%rates%actual_infection_rate, sim_stats%rates%actual_recovery_rate
            !$omp end critical(file_write)
         end if
      end if

   end subroutine run_simulation_loop

   subroutine close_output_files(stats_unit, events_unit, should_write_stats, should_write_events)
      integer(ik), intent(in), optional :: stats_unit, events_unit
      logical, intent(in) :: should_write_stats, should_write_events

      !$omp critical(file_write)
      if (should_write_stats) close(unit=stats_unit)
      if (should_write_events) close(unit=events_unit)
      !$omp end critical(file_write)
   end subroutine close_output_files

end program main
